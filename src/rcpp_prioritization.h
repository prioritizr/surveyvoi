#ifndef PRIORITIZATION_H
#define PRIORITIZATION_H

#include "package.h"
#include "functions.h"
#include "gurobi_c.h"

class Prioritization
{
public:
  // attributes
  std::size_t _n_f;     // number of features
  std::size_t _n_pu;    // number of planning units
  std::size_t _n_vars;  // number of variables in problem
  std::size_t _n_lcs;   // number of linear constraints in problem
  std::vector<int> _update_constraints_idx; // update indices for rij
  std::vector<int> _update_vars_idx;        // update indices for rij
  std::vector<double> _preweight_slope;  // objective function scaling parameter
  std::vector<double> _postweight_slope; // objective function scaling parameter
  GRBenv *_env;           // gurobi problem environment
  GRBmodel *_model;       // gurobi problem object

  // constructor
  Prioritization(std::size_t n_pu, std::size_t n_f,
                 Eigen::VectorXd &pu_costs,
                 Eigen::VectorXd &pu_locked_in,
                 Eigen::VectorXd &pu_locked_out,
                 int target, double budget, double gap) :
                 _env(NULL), _model(NULL) {
    // initialization
    _n_pu = n_pu;
    _n_f = n_f;
    _n_vars = n_pu + n_f + 1;
    _n_lcs = n_f + 1 + 1;

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
      lb[i] = -GRB_INFINITY;
    /// upper bounds
    std::vector<double> ub(_n_vars);
    for (std::size_t i = 0; i < _n_pu; ++i)
      ub[i] = 1.0 - pu_locked_out[i];
    for (std::size_t i = _n_pu; i < _n_vars; ++i)
      ub[i] = 0.0;

    /// linear component of objective function
    std::vector<double> obj(_n_vars);
    for (std::size_t i = 0; i < _n_pu; ++i)
      obj[i] = 0.0;
    for (std::size_t i = 0; i < _n_f; ++i)
      obj[_n_pu + i] = 0.0;
    obj[_n_pu + _n_f] = 1.0;

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
    GRBsetintparam(GRBgetenv(_model), "Presolve", -1); // enable presolve
    GRBsetintparam(GRBgetenv(_model), "Threads", 1);  // disable multi-threads
    GRBsetintparam(GRBgetenv(_model), "NumericFocus", 2);  // default precision
    GRBsetintparam(GRBgetenv(_model), "OutputFlag", 0);  // disable output
    GRBsetintparam(GRBgetenv(_model), "Seed", 1);  // disable output

    // model sense
    GRBsetintattr(_model, "ModelSense", GRB_MINIMIZE);

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
                    &feature_held_col_val[0], GRB_EQUAL, 0.0, NULL);
    }
    //// add constraints to calculate least represented feature
    //// initialize vectors
    std::vector<int> min_feature_held_col_idx(2);
    min_feature_held_col_idx[0] = 1.0;
    min_feature_held_col_idx[1] = _n_pu + _n_f;
    std::vector<double> min_feature_held_col_val(2);
    min_feature_held_col_val[0] = 1.0;
    min_feature_held_col_val[1] = -1.0;
    //// add constraints
    for (std::size_t j = 0; j < _n_f; ++j) {
      min_feature_held_col_idx[0] = _n_pu + j;
      GRBaddconstr(_model, 2, &min_feature_held_col_idx[0],
                    &min_feature_held_col_val[0], GRB_LESS_EQUAL, 0.0, NULL);
    }
    /// minimum number of planning units in solution constraint
    //// initialize vectors
    std::vector<int> min_pu_col_idx(_n_pu);
    std::iota(min_pu_col_idx.begin(), min_pu_col_idx.end(), 0);
    std::vector<double> min_pu_col_val(_n_pu, 1.0);
    //// add constraint
    GRBaddconstr(_model, _n_pu, &min_pu_col_idx[0],
                 min_pu_col_val.data(), GRB_GREATER_EQUAL, target, NULL);
    /// budget constraints
    //// add constraint
    GRBaddconstr(_model, _n_pu, &min_pu_col_idx[0],
                 pu_costs.data(), GRB_LESS_EQUAL, budget, NULL);
  };

  // destructor
  ~Prioritization() {
    GRBfreemodel(_model);
    GRBfreeenv(_env);
  };

  // add rij data to problem
  inline void add_rij_data(Eigen::MatrixXd &rij) {
    // update coefficients in the linear constraint matrix
    Eigen::MatrixXd rij2 = rij;
    rij2 = rij2.cwiseMin(1.0 - 1.0e-10);
    log_1m_matrix(rij2);
    rij2 = rij2.array() * 1.0e+3;
    GRBchgcoeffs(_model, _n_f * _n_pu, &_update_constraints_idx[0],
                 &_update_vars_idx[0], rij2.data());
    // return void
    return;
  }

  // solve problem
  inline void solve() {
    // generate solution
    GRBreset(_model, 0);
    int error = GRBoptimize(_model);
    // check for error
    if (static_cast<bool>(error)) {
      // throw error
      std::string err_message =
        "issue solving prioritization problem (code " +
        std::to_string(error) + ")";
      Rcpp::stop(err_message);
    }
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

  // extract solution
  inline void get_solution(Eigen::VectorXd &out) {
    // extract solution
    GRBgetdblattrarray(_model, GRB_DBL_ATTR_X, 0, _n_pu, out.data());
    // return void
    return;
  }

  // extract solution
  inline void get_solution(Eigen::ArrayXd &out) {
    // extract solution
    GRBgetdblattrarray(_model, GRB_DBL_ATTR_X, 0, _n_pu, out.data());
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
