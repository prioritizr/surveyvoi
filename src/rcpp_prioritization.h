#ifndef PRIORITIZATION_H
#define PRIORITIZATION_H

#include "package.h"
#include "gurobi_c.h"

class Prioritization
{
public:
  // attributes
  std::size_t _n_f;     // number of features
  std::size_t _n_pu;    // number of planning units
  std::size_t _n_vars;  // number of variables in problem
  std::size_t _n_lcs;   // number of linear constraints in problem
  std::size_t _n_approx_obj_fun_points; // number of approximation points
  double *_alpha; // objective function scaling parameter
  double *_gamma; // objective function scaling parameter
  std::vector<int> _update_constraints_idx; // update indices for rij
  std::vector<int> _update_vars_idx;        // update indices for rij

  GRBenv *_env;           // gurobi problem environment
  GRBmodel *_model;       // gurobi problem object

  // constructor
  Prioritization(std::size_t n_pu, std::size_t n_f,
                 Eigen::VectorXd &pu_costs,
                 Eigen::VectorXd &pu_locked_in,
                 Eigen::VectorXd &alpha,
                 Eigen::VectorXd &gamma,
                 std::size_t n_approx_obj_fun_points,
                 double budget, double gap) :
                 _env(NULL), _model(NULL) {
    // initialization
    _n_pu = n_pu;
    _n_f = n_f;
    _n_vars = n_pu + n_f;
    _n_lcs = n_f + 1;
    _n_approx_obj_fun_points = n_approx_obj_fun_points;
    _alpha = alpha.data();
    _gamma = gamma.data();

    // create update indices for rij data
    /// add constraint matrix row numbers
    _update_constraints_idx.reserve(_n_f * _n_pu);
    for (std::size_t i = 0; i < _n_pu; ++i)
      for (std::size_t j = 0; j < _n_f; ++j)
        _update_constraints_idx.push_back(j);
    /// add constraint matrix column numbers
    _update_vars_idx.reserve(_n_f * _n_pu);
    for (std::size_t i = 0; i < _n_pu; ++i)
      for (std::size_t j = 0; j < _n_f; ++j)
        _update_vars_idx.push_back(i);

    // model variables
    /// lower bounds
    std::vector<double> lb(_n_vars);
    for (std::size_t i = 0; i < _n_pu; ++i)
      lb[i] = pu_locked_in[i];
    for (std::size_t i = _n_pu; i < _n_vars; ++i)
      lb[i] = 0.0;
    /// upper bounds
    std::vector<double> ub(_n_vars);
    for (std::size_t i = 0; i < _n_pu; ++i)
      ub[i] = 1.0;
    for (std::size_t i = _n_pu; i < _n_vars; ++i)
      ub[i] = GRB_INFINITY;
    /// linear component of objective function
    std::vector<double> obj(_n_vars, 0.0);
    /// variable types
    std::vector<char> vtype(_n_vars);
    for (std::size_t i = 0; i < _n_pu; ++i)
      vtype[i] = 'B';
    for (std::size_t i = _n_pu; i < _n_vars; ++i)
      vtype[i] = 'C';

    // initialize environmment
    GRBloadenv(&_env, NULL);

    // initialize model
    GRBnewmodel(_env, &_model, "v", _n_vars, &obj[0], &lb[0], &ub[0],
                &vtype[0], NULL);

    // set parameters for solving problem
    GRBsetdblparam(GRBgetenv(_model), "MIPGap", gap); // optimal solutions
    GRBsetintparam(GRBgetenv(_model), "Presolve", 2); // enable presolve
    GRBsetintparam(GRBgetenv(_model), "Threads", 1);  // disable multi-threads
    GRBsetintparam(GRBgetenv(_model), "NumericFocus", 0);  // default precision
    GRBsetintparam(GRBgetenv(_model), "OutputFlag", 0);  // disable output

    // model sense
    GRBsetintattr(_model, "ModelSense", GRB_MAXIMIZE);

    // add constraints to problem
    /// feature held calculation constraints
    //// initialize vectors
    std::vector<int> feature_held_col_idx(_n_pu + 1);
    std::iota(feature_held_col_idx.begin(), feature_held_col_idx.end(), 0);
    std::vector<double> feature_held_col_val(_n_pu + 1, 1.0);
    feature_held_col_val[_n_pu] = -1.0;
    //// add constraints
    for (std::size_t j = 0; j < _n_f; ++j) {
      feature_held_col_idx[_n_pu] = _n_pu + j;
      GRBaddconstr(_model, _n_pu + 1, &feature_held_col_idx[0],
                    &feature_held_col_val[0], '=', 0.0, NULL);
    }
    /// budget constraints
    //// initialize vectors
    std::vector<int> budget_col_idx(_n_pu);
    std::iota(budget_col_idx.begin(), budget_col_idx.end(), 0);
    //// add constraint
    GRBaddconstr(_model, _n_pu, &feature_held_col_idx[0],
                 pu_costs.data(), '<', budget, NULL);
  };

  // destructor
  ~Prioritization() {
    GRBfreemodel(_model);
    GRBfreeenv(_env);

  };

  // add rij data to problem
  inline void add_rij_data(Eigen::MatrixXd &rij) {
    // update coefficients in the linear constraint matrix
    GRBchgcoeffs(_model, _n_f * _n_pu, &_update_constraints_idx[0],
                 &_update_vars_idx[0], rij.data());
    // add piece-wise linear objective function components
    /// initialize dummmy variables if a feature only contains zeros
    Eigen::VectorXd dummy_held(3);
    Eigen::VectorXd dummy_benefit(3);
    dummy_held[0] = 0.0;
    dummy_benefit[0] = 0.0;
    dummy_held[1] = 0.5;
    dummy_benefit[1] = 0.5;
    /// initialize variables
    double curr_value;
    Eigen::VectorXd min_feature_held = rij.rowwise().minCoeff();
    min_feature_held *= 0.99;
    Eigen::VectorXd sum_feature_held = rij.rowwise().sum();
    sum_feature_held *= 1.01;
    Eigen::VectorXd incr_feature_held =
      (sum_feature_held.array() - min_feature_held.array()) /
      static_cast<double>(_n_approx_obj_fun_points - 2);
    Eigen::VectorXd obj_feature_held(_n_approx_obj_fun_points);
    Eigen::VectorXd obj_feature_benefit(_n_approx_obj_fun_points);
    obj_feature_held[0] = 0.0;
    obj_feature_benefit[0] = 0.0;
    /// add piece-wise components
    for (std::size_t j = 0; j < _n_f; ++j) {
      /// calculate held values
      curr_value = min_feature_held(j);
      for (std::size_t k = 1; k < _n_approx_obj_fun_points; ++k) {
        obj_feature_held[k] = curr_value;
        curr_value += incr_feature_held[j];
      }
      /// calculate benefit values
      obj_feature_benefit =
        (*(_alpha + j) * obj_feature_held.array()).pow(*(_gamma + j));
      obj_feature_benefit[0] = 0.0;
      /// add component
      if ((obj_feature_benefit.maxCoeff() - obj_feature_benefit[0]) > 1.0e-5) {
        // use actual peice-wise linear objectives if the feature doesn't
        // have zeros in all planning units
        GRBsetpwlobj(_model, _n_pu + j, _n_approx_obj_fun_points,
                     obj_feature_held.data(),
                     obj_feature_benefit.data());
      } else {
        // use dummy peice-wise linear objectives if there is zero variation
        // in the benefit function
        if (obj_feature_benefit[1] > 1.0e-5) {
          // if the benefit function has constant values that are greater
          // than zero, then set this value as the maximum value
          dummy_held[2] = obj_feature_held[1];
          dummy_benefit[2] = obj_feature_benefit[1];
        } else {
          dummy_held[2] = 1.0;
          dummy_benefit[2] = 1.0;
        }
        GRBsetpwlobj(_model, _n_pu + j, 3, dummy_held.data(),
                     dummy_benefit.data());
      }
    }
    // return void
    return;
  }

  // solve problem
  inline void solve() {
    // generate solution
    int error = GRBoptimize(_model);
    // check for error
    if (static_cast<bool>(error))
      Rcpp::stop("issue solving prioritization problem");
    // return void
    return;
  }

  // extract solution
  inline void get_solution(std::vector<bool> &out) {
    // extract solution
    Eigen::VectorXd solution(_n_pu);
    GRBgetdblattrarray(_model, GRB_DBL_ATTR_X, 0, _n_pu, solution.data());
    // store it in the output
    for (std::size_t i = 0; i < _n_pu; ++i)
      out[i] = solution[i] > 0.5;
    // return void
    return;
  }

  // objective value
  inline double get_obj_value() {
    double out;
    GRBgetdblattr(_model, GRB_DBL_ATTR_OBJVAL, &out);
    return out;
  }

  // save problem to file
  inline void save_problem(std::string &file_path) {
    GRBwrite(_model, file_path.c_str());
  }

};

#endif
