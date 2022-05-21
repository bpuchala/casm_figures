#include <cstddef>
#include "casm/clexulator/BaseClexulator.hh"
#include "casm/clexulator/BasicClexParamPack.hh"
#include "casm/global/eigen.hh"
#include "casm/clexulator/BasicClexParamPack.hh"




/****** PROJECT SPECIFICATIONS ******

         ****** prim.json ******

{
  "basis" : [
    {
      "coordinate" : [ 0.000000000000, 0.000000000000, 0.000000000000 ],
      "occupants" : [ "A", "B" ]
    }
  ],
  "coordinate_mode" : "Fractional",
  "lattice_vectors" : [
    [ 0.000000000000, 0.500000000000, 0.500000000000 ],
    [ 0.500000000000, 0.000000000000, 0.500000000000 ],
    [ 0.500000000000, 0.500000000000, 0.000000000000 ]
  ],
  "title" : "fcc_binary"
}

        ****** bspecs.json ******

{
  "basis_function_specs" : {
    "dof_specs" : {
      "occ" : {
        "site_basis_functions" : "CHEBYCHEV"
      }
    },
    "dofs" : [ "occ" ],
    "global_max_poly_order" : -1,
    "param_pack_type" : "DEFAULT"
  },
  "cluster_specs" : {
    "method" : "periodic_max_length",
    "params" : {
      "generating_group" : [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47 ],
      "orbit_branch_specs" : {
        "0" : {
          "max_length" : 0.000000000000
        },
        "1" : {
          "max_length" : 0.000000000000
        },
        "2" : {
          "max_length" : 3.010000000000
        },
        "3" : {
          "max_length" : 1.510000000000
        },
        "4" : {
          "max_length" : 1.010000000000
        }
      }
    }
  }
}

**/


/// \brief Returns a clexulator::BaseClexulator* owning a fcc_binary_Clexulator_default
extern "C" CASM::clexulator::BaseClexulator *make_fcc_binary_Clexulator_default();

namespace CASM {
namespace clexulator {

  /****** GENERATED CLEXPARAMPACK DEFINITION ******/


  typedef BasicClexParamPack ParamPack;


  /****** GENERATED CLEXULATOR DEFINITION ******/

  class fcc_binary_Clexulator_default : public clexulator::BaseClexulator {

  public:

    fcc_binary_Clexulator_default();

    ~fcc_binary_Clexulator_default();

    ClexParamPack const &param_pack() const override {
      return m_params;
    }

    ClexParamPack &param_pack() override {
      return m_params;
    }


    template<typename Scalar>
    Scalar eval_bfunc_0_0() const;

    template<typename Scalar>
    Scalar eval_bfunc_1_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_1_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_1_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_2_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_2_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_2_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_3_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_3_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_3_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_4_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_4_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_4_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_5_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_5_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_5_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_6_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_6_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_6_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_7_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_7_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_7_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_8_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_8_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_8_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_9_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_9_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_9_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_10_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_10_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_10_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_11_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_11_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_11_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_12_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_12_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_12_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_13_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_13_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_13_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_14_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_14_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_14_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_15_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_15_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_15_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_16_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_16_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_16_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_17_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_17_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_17_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_18_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_18_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_18_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_19_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_19_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_19_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_20_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_20_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_20_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_21_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_21_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_21_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_22_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_22_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_22_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_23_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_23_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_23_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_24_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_24_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_24_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_25_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_25_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_25_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_26_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_26_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_26_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_27_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_27_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_27_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_28_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_28_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_28_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_29_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_29_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_29_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_30_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_30_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_30_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_31_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_31_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_31_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_32_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_32_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_32_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_33_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_33_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_33_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_34_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_34_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_34_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_35_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_35_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_35_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_36_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_36_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_36_0_at_0(int occ_i, int occ_f) const;

    template<typename Scalar>
    Scalar eval_bfunc_37_0() const;

    template<typename Scalar>
    Scalar site_eval_bfunc_37_0_at_0() const;

    template<typename Scalar>
    Scalar site_deval_bfunc_37_0_at_0(int occ_i, int occ_f) const;


  private:

    // ParamPack object, which stores temporary data for calculations
    mutable ParamPack m_params;

    // typedef for method pointers of scalar type double
    typedef double (fcc_binary_Clexulator_default::*BasisFuncPtr_0)() const;

    // typedef for method pointers
    typedef double (fcc_binary_Clexulator_default::*DeltaBasisFuncPtr_0)(int, int) const;

    // array of pointers to member functions for calculating basis functions of scalar type double
    BasisFuncPtr_0 m_orbit_func_table_0[38];

    // array of pointers to member functions for calculating flower functions of scalar type double
    BasisFuncPtr_0 m_flower_func_table_0[1][38];

    // array of pointers to member functions for calculating DELTA flower functions of scalar type double
    DeltaBasisFuncPtr_0 m_delta_func_table_0[1][38];

    // Occupation Function tables for basis sites in asymmetric unit 0:
    //   - basis site 0:
    double m_occ_func_0_0[2];

    //ClexParamPack allocation for evaluated correlations 
    ParamPack::Key m_corr_param_key;
    //ClexParamPack allocation for DoF occ
    ParamPack::Key m_occ_site_func_param_key;

    /// \brief Clone the fcc_binary_Clexulator_default
    BaseClexulator *_clone() const override {
      return new fcc_binary_Clexulator_default(*this);
    }

    /// \brief Calculate contribution to global correlations from one unit cell
    /// Result is recorded in ClexParamPack
    void _calc_global_corr_contribution() const override;

    /// \brief Calculate contribution to global correlations from one unit cell     /// Result is recorded in double array starting at corr_begin
    void _calc_global_corr_contribution(double *corr_begin) const override;

    /// \brief Calculate contribution to select global correlations from one unit cell into ClexParamPack
    /// Result is recorded in ClexParamPack
    void _calc_restricted_global_corr_contribution(size_type const *ind_list_begin, size_type const *ind_list_end) const override;

    /// \brief Calculate contribution to select global correlations from one unit cell
    /// Result is recorded in double array starting at corr_begin
    void _calc_restricted_global_corr_contribution(double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const override;

    /// \brief Calculate point correlations about neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in ClexParamPack
    void _calc_point_corr(int nlist_ind) const override;

    /// \brief Calculate point correlations about neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in double array starting at corr_begin
    void _calc_point_corr(int nlist_ind, double *corr_begin) const override;

    /// \brief Calculate select point correlations about neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in ClexParamPack
    void _calc_restricted_point_corr(int nlist_ind, size_type const *ind_list_begin, size_type const *ind_list_end) const override;

    /// \brief Calculate select point correlations about neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in double array starting at corr_begin
    void _calc_restricted_point_corr(int nlist_ind, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const override;

    /// \brief Calculate the change in point correlations due to changing an occupant at neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in ClexParamPack
    void _calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f) const override;

    /// \brief Calculate the change in point correlations due to changing an occupant at neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in double array starting at corr_begin
    void _calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double *corr_begin) const override;

    /// \brief Calculate the change in select point correlations due to changing an occupant at neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in ClexParamPack
    void _calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f, size_type const *ind_list_begin, size_type const *ind_list_end) const override;

    /// \brief Calculate the change in select point correlations due to changing an occupant at neighbor site 'nlist_ind'
    /// For global clexulators, 'nlist_ind' only ranges over sites in the cell
    /// For local clexulators, 'nlist_ind' ranges over all sites in the neighborhood
    /// Result is recorded in double array starting at corr_begin
    void _calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const override;

    template<typename Scalar>
    void _global_prepare() const;

    template<typename Scalar>
    void _point_prepare(int nlist_ind) const;

    // Occupation Function evaluators and accessors for basis site 0:
    double const &eval_occ_func_0_0(const int &nlist_ind) const {
      return m_occ_func_0_0[_occ(nlist_ind)];
    }

    double const &occ_func_0_0(const int &nlist_ind) const {
      return m_params.read(m_occ_site_func_param_key, 0, nlist_ind);
    }

    //default functions for basis function evaluation
    template <typename Scalar>
    Scalar zero_func() const {
      return Scalar(0.0);
    }

    template <typename Scalar>
    Scalar zero_func(int, int) const {
      return Scalar(0.0);
    }


  };

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  fcc_binary_Clexulator_default::fcc_binary_Clexulator_default() :
    BaseClexulator(459, 38, 1) {
    m_occ_func_0_0[0] = -1.0000000000, m_occ_func_0_0[1] = 1.0000000000;

    m_occ_site_func_param_key = m_params.allocate("occ_site_func", 1, 459, true);

    m_corr_param_key = m_params.allocate("corr", corr_size(), 1, false);

    m_orbit_func_table_0[0] = &fcc_binary_Clexulator_default::eval_bfunc_0_0<double>;
    m_orbit_func_table_0[1] = &fcc_binary_Clexulator_default::eval_bfunc_1_0<double>;
    m_orbit_func_table_0[2] = &fcc_binary_Clexulator_default::eval_bfunc_2_0<double>;
    m_orbit_func_table_0[3] = &fcc_binary_Clexulator_default::eval_bfunc_3_0<double>;
    m_orbit_func_table_0[4] = &fcc_binary_Clexulator_default::eval_bfunc_4_0<double>;
    m_orbit_func_table_0[5] = &fcc_binary_Clexulator_default::eval_bfunc_5_0<double>;
    m_orbit_func_table_0[6] = &fcc_binary_Clexulator_default::eval_bfunc_6_0<double>;
    m_orbit_func_table_0[7] = &fcc_binary_Clexulator_default::eval_bfunc_7_0<double>;
    m_orbit_func_table_0[8] = &fcc_binary_Clexulator_default::eval_bfunc_8_0<double>;
    m_orbit_func_table_0[9] = &fcc_binary_Clexulator_default::eval_bfunc_9_0<double>;
    m_orbit_func_table_0[10] = &fcc_binary_Clexulator_default::eval_bfunc_10_0<double>;
    m_orbit_func_table_0[11] = &fcc_binary_Clexulator_default::eval_bfunc_11_0<double>;
    m_orbit_func_table_0[12] = &fcc_binary_Clexulator_default::eval_bfunc_12_0<double>;
    m_orbit_func_table_0[13] = &fcc_binary_Clexulator_default::eval_bfunc_13_0<double>;
    m_orbit_func_table_0[14] = &fcc_binary_Clexulator_default::eval_bfunc_14_0<double>;
    m_orbit_func_table_0[15] = &fcc_binary_Clexulator_default::eval_bfunc_15_0<double>;
    m_orbit_func_table_0[16] = &fcc_binary_Clexulator_default::eval_bfunc_16_0<double>;
    m_orbit_func_table_0[17] = &fcc_binary_Clexulator_default::eval_bfunc_17_0<double>;
    m_orbit_func_table_0[18] = &fcc_binary_Clexulator_default::eval_bfunc_18_0<double>;
    m_orbit_func_table_0[19] = &fcc_binary_Clexulator_default::eval_bfunc_19_0<double>;
    m_orbit_func_table_0[20] = &fcc_binary_Clexulator_default::eval_bfunc_20_0<double>;
    m_orbit_func_table_0[21] = &fcc_binary_Clexulator_default::eval_bfunc_21_0<double>;
    m_orbit_func_table_0[22] = &fcc_binary_Clexulator_default::eval_bfunc_22_0<double>;
    m_orbit_func_table_0[23] = &fcc_binary_Clexulator_default::eval_bfunc_23_0<double>;
    m_orbit_func_table_0[24] = &fcc_binary_Clexulator_default::eval_bfunc_24_0<double>;
    m_orbit_func_table_0[25] = &fcc_binary_Clexulator_default::eval_bfunc_25_0<double>;
    m_orbit_func_table_0[26] = &fcc_binary_Clexulator_default::eval_bfunc_26_0<double>;
    m_orbit_func_table_0[27] = &fcc_binary_Clexulator_default::eval_bfunc_27_0<double>;
    m_orbit_func_table_0[28] = &fcc_binary_Clexulator_default::eval_bfunc_28_0<double>;
    m_orbit_func_table_0[29] = &fcc_binary_Clexulator_default::eval_bfunc_29_0<double>;
    m_orbit_func_table_0[30] = &fcc_binary_Clexulator_default::eval_bfunc_30_0<double>;
    m_orbit_func_table_0[31] = &fcc_binary_Clexulator_default::eval_bfunc_31_0<double>;
    m_orbit_func_table_0[32] = &fcc_binary_Clexulator_default::eval_bfunc_32_0<double>;
    m_orbit_func_table_0[33] = &fcc_binary_Clexulator_default::eval_bfunc_33_0<double>;
    m_orbit_func_table_0[34] = &fcc_binary_Clexulator_default::eval_bfunc_34_0<double>;
    m_orbit_func_table_0[35] = &fcc_binary_Clexulator_default::eval_bfunc_35_0<double>;
    m_orbit_func_table_0[36] = &fcc_binary_Clexulator_default::eval_bfunc_36_0<double>;
    m_orbit_func_table_0[37] = &fcc_binary_Clexulator_default::eval_bfunc_37_0<double>;


    m_flower_func_table_0[0][0] = &fcc_binary_Clexulator_default::zero_func<double>;
    m_flower_func_table_0[0][1] = &fcc_binary_Clexulator_default::site_eval_bfunc_1_0_at_0<double>;
    m_flower_func_table_0[0][2] = &fcc_binary_Clexulator_default::site_eval_bfunc_2_0_at_0<double>;
    m_flower_func_table_0[0][3] = &fcc_binary_Clexulator_default::site_eval_bfunc_3_0_at_0<double>;
    m_flower_func_table_0[0][4] = &fcc_binary_Clexulator_default::site_eval_bfunc_4_0_at_0<double>;
    m_flower_func_table_0[0][5] = &fcc_binary_Clexulator_default::site_eval_bfunc_5_0_at_0<double>;
    m_flower_func_table_0[0][6] = &fcc_binary_Clexulator_default::site_eval_bfunc_6_0_at_0<double>;
    m_flower_func_table_0[0][7] = &fcc_binary_Clexulator_default::site_eval_bfunc_7_0_at_0<double>;
    m_flower_func_table_0[0][8] = &fcc_binary_Clexulator_default::site_eval_bfunc_8_0_at_0<double>;
    m_flower_func_table_0[0][9] = &fcc_binary_Clexulator_default::site_eval_bfunc_9_0_at_0<double>;
    m_flower_func_table_0[0][10] = &fcc_binary_Clexulator_default::site_eval_bfunc_10_0_at_0<double>;
    m_flower_func_table_0[0][11] = &fcc_binary_Clexulator_default::site_eval_bfunc_11_0_at_0<double>;
    m_flower_func_table_0[0][12] = &fcc_binary_Clexulator_default::site_eval_bfunc_12_0_at_0<double>;
    m_flower_func_table_0[0][13] = &fcc_binary_Clexulator_default::site_eval_bfunc_13_0_at_0<double>;
    m_flower_func_table_0[0][14] = &fcc_binary_Clexulator_default::site_eval_bfunc_14_0_at_0<double>;
    m_flower_func_table_0[0][15] = &fcc_binary_Clexulator_default::site_eval_bfunc_15_0_at_0<double>;
    m_flower_func_table_0[0][16] = &fcc_binary_Clexulator_default::site_eval_bfunc_16_0_at_0<double>;
    m_flower_func_table_0[0][17] = &fcc_binary_Clexulator_default::site_eval_bfunc_17_0_at_0<double>;
    m_flower_func_table_0[0][18] = &fcc_binary_Clexulator_default::site_eval_bfunc_18_0_at_0<double>;
    m_flower_func_table_0[0][19] = &fcc_binary_Clexulator_default::site_eval_bfunc_19_0_at_0<double>;
    m_flower_func_table_0[0][20] = &fcc_binary_Clexulator_default::site_eval_bfunc_20_0_at_0<double>;
    m_flower_func_table_0[0][21] = &fcc_binary_Clexulator_default::site_eval_bfunc_21_0_at_0<double>;
    m_flower_func_table_0[0][22] = &fcc_binary_Clexulator_default::site_eval_bfunc_22_0_at_0<double>;
    m_flower_func_table_0[0][23] = &fcc_binary_Clexulator_default::site_eval_bfunc_23_0_at_0<double>;
    m_flower_func_table_0[0][24] = &fcc_binary_Clexulator_default::site_eval_bfunc_24_0_at_0<double>;
    m_flower_func_table_0[0][25] = &fcc_binary_Clexulator_default::site_eval_bfunc_25_0_at_0<double>;
    m_flower_func_table_0[0][26] = &fcc_binary_Clexulator_default::site_eval_bfunc_26_0_at_0<double>;
    m_flower_func_table_0[0][27] = &fcc_binary_Clexulator_default::site_eval_bfunc_27_0_at_0<double>;
    m_flower_func_table_0[0][28] = &fcc_binary_Clexulator_default::site_eval_bfunc_28_0_at_0<double>;
    m_flower_func_table_0[0][29] = &fcc_binary_Clexulator_default::site_eval_bfunc_29_0_at_0<double>;
    m_flower_func_table_0[0][30] = &fcc_binary_Clexulator_default::site_eval_bfunc_30_0_at_0<double>;
    m_flower_func_table_0[0][31] = &fcc_binary_Clexulator_default::site_eval_bfunc_31_0_at_0<double>;
    m_flower_func_table_0[0][32] = &fcc_binary_Clexulator_default::site_eval_bfunc_32_0_at_0<double>;
    m_flower_func_table_0[0][33] = &fcc_binary_Clexulator_default::site_eval_bfunc_33_0_at_0<double>;
    m_flower_func_table_0[0][34] = &fcc_binary_Clexulator_default::site_eval_bfunc_34_0_at_0<double>;
    m_flower_func_table_0[0][35] = &fcc_binary_Clexulator_default::site_eval_bfunc_35_0_at_0<double>;
    m_flower_func_table_0[0][36] = &fcc_binary_Clexulator_default::site_eval_bfunc_36_0_at_0<double>;
    m_flower_func_table_0[0][37] = &fcc_binary_Clexulator_default::site_eval_bfunc_37_0_at_0<double>;


    m_delta_func_table_0[0][0] = &fcc_binary_Clexulator_default::zero_func<double>;
    m_delta_func_table_0[0][1] = &fcc_binary_Clexulator_default::site_deval_bfunc_1_0_at_0<double>;
    m_delta_func_table_0[0][2] = &fcc_binary_Clexulator_default::site_deval_bfunc_2_0_at_0<double>;
    m_delta_func_table_0[0][3] = &fcc_binary_Clexulator_default::site_deval_bfunc_3_0_at_0<double>;
    m_delta_func_table_0[0][4] = &fcc_binary_Clexulator_default::site_deval_bfunc_4_0_at_0<double>;
    m_delta_func_table_0[0][5] = &fcc_binary_Clexulator_default::site_deval_bfunc_5_0_at_0<double>;
    m_delta_func_table_0[0][6] = &fcc_binary_Clexulator_default::site_deval_bfunc_6_0_at_0<double>;
    m_delta_func_table_0[0][7] = &fcc_binary_Clexulator_default::site_deval_bfunc_7_0_at_0<double>;
    m_delta_func_table_0[0][8] = &fcc_binary_Clexulator_default::site_deval_bfunc_8_0_at_0<double>;
    m_delta_func_table_0[0][9] = &fcc_binary_Clexulator_default::site_deval_bfunc_9_0_at_0<double>;
    m_delta_func_table_0[0][10] = &fcc_binary_Clexulator_default::site_deval_bfunc_10_0_at_0<double>;
    m_delta_func_table_0[0][11] = &fcc_binary_Clexulator_default::site_deval_bfunc_11_0_at_0<double>;
    m_delta_func_table_0[0][12] = &fcc_binary_Clexulator_default::site_deval_bfunc_12_0_at_0<double>;
    m_delta_func_table_0[0][13] = &fcc_binary_Clexulator_default::site_deval_bfunc_13_0_at_0<double>;
    m_delta_func_table_0[0][14] = &fcc_binary_Clexulator_default::site_deval_bfunc_14_0_at_0<double>;
    m_delta_func_table_0[0][15] = &fcc_binary_Clexulator_default::site_deval_bfunc_15_0_at_0<double>;
    m_delta_func_table_0[0][16] = &fcc_binary_Clexulator_default::site_deval_bfunc_16_0_at_0<double>;
    m_delta_func_table_0[0][17] = &fcc_binary_Clexulator_default::site_deval_bfunc_17_0_at_0<double>;
    m_delta_func_table_0[0][18] = &fcc_binary_Clexulator_default::site_deval_bfunc_18_0_at_0<double>;
    m_delta_func_table_0[0][19] = &fcc_binary_Clexulator_default::site_deval_bfunc_19_0_at_0<double>;
    m_delta_func_table_0[0][20] = &fcc_binary_Clexulator_default::site_deval_bfunc_20_0_at_0<double>;
    m_delta_func_table_0[0][21] = &fcc_binary_Clexulator_default::site_deval_bfunc_21_0_at_0<double>;
    m_delta_func_table_0[0][22] = &fcc_binary_Clexulator_default::site_deval_bfunc_22_0_at_0<double>;
    m_delta_func_table_0[0][23] = &fcc_binary_Clexulator_default::site_deval_bfunc_23_0_at_0<double>;
    m_delta_func_table_0[0][24] = &fcc_binary_Clexulator_default::site_deval_bfunc_24_0_at_0<double>;
    m_delta_func_table_0[0][25] = &fcc_binary_Clexulator_default::site_deval_bfunc_25_0_at_0<double>;
    m_delta_func_table_0[0][26] = &fcc_binary_Clexulator_default::site_deval_bfunc_26_0_at_0<double>;
    m_delta_func_table_0[0][27] = &fcc_binary_Clexulator_default::site_deval_bfunc_27_0_at_0<double>;
    m_delta_func_table_0[0][28] = &fcc_binary_Clexulator_default::site_deval_bfunc_28_0_at_0<double>;
    m_delta_func_table_0[0][29] = &fcc_binary_Clexulator_default::site_deval_bfunc_29_0_at_0<double>;
    m_delta_func_table_0[0][30] = &fcc_binary_Clexulator_default::site_deval_bfunc_30_0_at_0<double>;
    m_delta_func_table_0[0][31] = &fcc_binary_Clexulator_default::site_deval_bfunc_31_0_at_0<double>;
    m_delta_func_table_0[0][32] = &fcc_binary_Clexulator_default::site_deval_bfunc_32_0_at_0<double>;
    m_delta_func_table_0[0][33] = &fcc_binary_Clexulator_default::site_deval_bfunc_33_0_at_0<double>;
    m_delta_func_table_0[0][34] = &fcc_binary_Clexulator_default::site_deval_bfunc_34_0_at_0<double>;
    m_delta_func_table_0[0][35] = &fcc_binary_Clexulator_default::site_deval_bfunc_35_0_at_0<double>;
    m_delta_func_table_0[0][36] = &fcc_binary_Clexulator_default::site_deval_bfunc_36_0_at_0<double>;
    m_delta_func_table_0[0][37] = &fcc_binary_Clexulator_default::site_deval_bfunc_37_0_at_0<double>;


    m_weight_matrix.row(0) << 2, 1, 1;
    m_weight_matrix.row(1) << 1, 2, 1;
    m_weight_matrix.row(2) << 1, 1, 2;

    m_sublat_indices = std::set<int>{0};

    m_n_sublattices = 1;

    m_neighborhood = std::set<xtal::UnitCell> {
      xtal::UnitCell(-5, 1, 1),
      xtal::UnitCell(-5, 1, 2),
      xtal::UnitCell(-5, 1, 3),
      xtal::UnitCell(-5, 2, 1),
      xtal::UnitCell(-5, 2, 2),
      xtal::UnitCell(-5, 3, 1),
      xtal::UnitCell(-4, -1, 1),
      xtal::UnitCell(-4, -1, 2),
      xtal::UnitCell(-4, -1, 3),
      xtal::UnitCell(-4, -1, 4),
      xtal::UnitCell(-4, 0, 0),
      xtal::UnitCell(-4, 0, 1),
      xtal::UnitCell(-4, 0, 2),
      xtal::UnitCell(-4, 0, 3),
      xtal::UnitCell(-4, 0, 4),
      xtal::UnitCell(-4, 1, -1),
      xtal::UnitCell(-4, 1, 0),
      xtal::UnitCell(-4, 1, 1),
      xtal::UnitCell(-4, 1, 2),
      xtal::UnitCell(-4, 1, 3),
      xtal::UnitCell(-4, 1, 4),
      xtal::UnitCell(-4, 2, -1),
      xtal::UnitCell(-4, 2, 0),
      xtal::UnitCell(-4, 2, 1),
      xtal::UnitCell(-4, 2, 2),
      xtal::UnitCell(-4, 2, 3),
      xtal::UnitCell(-4, 3, -1),
      xtal::UnitCell(-4, 3, 0),
      xtal::UnitCell(-4, 3, 1),
      xtal::UnitCell(-4, 3, 2),
      xtal::UnitCell(-4, 4, -1),
      xtal::UnitCell(-4, 4, 0),
      xtal::UnitCell(-4, 4, 1),
      xtal::UnitCell(-3, -3, 3),
      xtal::UnitCell(-3, -2, 1),
      xtal::UnitCell(-3, -2, 2),
      xtal::UnitCell(-3, -2, 3),
      xtal::UnitCell(-3, -2, 4),
      xtal::UnitCell(-3, -1, -1),
      xtal::UnitCell(-3, -1, 0),
      xtal::UnitCell(-3, -1, 1),
      xtal::UnitCell(-3, -1, 2),
      xtal::UnitCell(-3, -1, 3),
      xtal::UnitCell(-3, -1, 4),
      xtal::UnitCell(-3, -1, 5),
      xtal::UnitCell(-3, 0, -1),
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(-3, 0, 1),
      xtal::UnitCell(-3, 0, 2),
      xtal::UnitCell(-3, 0, 3),
      xtal::UnitCell(-3, 0, 4),
      xtal::UnitCell(-3, 1, -2),
      xtal::UnitCell(-3, 1, -1),
      xtal::UnitCell(-3, 1, 0),
      xtal::UnitCell(-3, 1, 1),
      xtal::UnitCell(-3, 1, 2),
      xtal::UnitCell(-3, 1, 3),
      xtal::UnitCell(-3, 1, 4),
      xtal::UnitCell(-3, 2, -2),
      xtal::UnitCell(-3, 2, -1),
      xtal::UnitCell(-3, 2, 0),
      xtal::UnitCell(-3, 2, 1),
      xtal::UnitCell(-3, 2, 2),
      xtal::UnitCell(-3, 2, 3),
      xtal::UnitCell(-3, 3, -3),
      xtal::UnitCell(-3, 3, -2),
      xtal::UnitCell(-3, 3, -1),
      xtal::UnitCell(-3, 3, 0),
      xtal::UnitCell(-3, 3, 1),
      xtal::UnitCell(-3, 3, 2),
      xtal::UnitCell(-3, 3, 3),
      xtal::UnitCell(-3, 4, -2),
      xtal::UnitCell(-3, 4, -1),
      xtal::UnitCell(-3, 4, 0),
      xtal::UnitCell(-3, 4, 1),
      xtal::UnitCell(-3, 5, -1),
      xtal::UnitCell(-2, -3, 1),
      xtal::UnitCell(-2, -3, 2),
      xtal::UnitCell(-2, -3, 3),
      xtal::UnitCell(-2, -3, 4),
      xtal::UnitCell(-2, -2, -1),
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, -2, 1),
      xtal::UnitCell(-2, -2, 2),
      xtal::UnitCell(-2, -2, 3),
      xtal::UnitCell(-2, -2, 4),
      xtal::UnitCell(-2, -2, 5),
      xtal::UnitCell(-2, -1, -2),
      xtal::UnitCell(-2, -1, -1),
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, -1, 1),
      xtal::UnitCell(-2, -1, 2),
      xtal::UnitCell(-2, -1, 3),
      xtal::UnitCell(-2, -1, 4),
      xtal::UnitCell(-2, -1, 5),
      xtal::UnitCell(-2, 0, -2),
      xtal::UnitCell(-2, 0, -1),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 0, 3),
      xtal::UnitCell(-2, 0, 4),
      xtal::UnitCell(-2, 1, -3),
      xtal::UnitCell(-2, 1, -2),
      xtal::UnitCell(-2, 1, -1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-2, 1, 2),
      xtal::UnitCell(-2, 1, 3),
      xtal::UnitCell(-2, 1, 4),
      xtal::UnitCell(-2, 2, -3),
      xtal::UnitCell(-2, 2, -2),
      xtal::UnitCell(-2, 2, -1),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-2, 2, 1),
      xtal::UnitCell(-2, 2, 2),
      xtal::UnitCell(-2, 2, 3),
      xtal::UnitCell(-2, 3, -3),
      xtal::UnitCell(-2, 3, -2),
      xtal::UnitCell(-2, 3, -1),
      xtal::UnitCell(-2, 3, 0),
      xtal::UnitCell(-2, 3, 1),
      xtal::UnitCell(-2, 3, 2),
      xtal::UnitCell(-2, 4, -3),
      xtal::UnitCell(-2, 4, -2),
      xtal::UnitCell(-2, 4, -1),
      xtal::UnitCell(-2, 4, 0),
      xtal::UnitCell(-2, 4, 1),
      xtal::UnitCell(-2, 5, -2),
      xtal::UnitCell(-2, 5, -1),
      xtal::UnitCell(-1, -4, 1),
      xtal::UnitCell(-1, -4, 2),
      xtal::UnitCell(-1, -4, 3),
      xtal::UnitCell(-1, -4, 4),
      xtal::UnitCell(-1, -3, -1),
      xtal::UnitCell(-1, -3, 0),
      xtal::UnitCell(-1, -3, 1),
      xtal::UnitCell(-1, -3, 2),
      xtal::UnitCell(-1, -3, 3),
      xtal::UnitCell(-1, -3, 4),
      xtal::UnitCell(-1, -3, 5),
      xtal::UnitCell(-1, -2, -2),
      xtal::UnitCell(-1, -2, -1),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -2, 1),
      xtal::UnitCell(-1, -2, 2),
      xtal::UnitCell(-1, -2, 3),
      xtal::UnitCell(-1, -2, 4),
      xtal::UnitCell(-1, -2, 5),
      xtal::UnitCell(-1, -1, -3),
      xtal::UnitCell(-1, -1, -2),
      xtal::UnitCell(-1, -1, -1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, -1, 3),
      xtal::UnitCell(-1, -1, 4),
      xtal::UnitCell(-1, -1, 5),
      xtal::UnitCell(-1, 0, -3),
      xtal::UnitCell(-1, 0, -2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 0, 3),
      xtal::UnitCell(-1, 0, 4),
      xtal::UnitCell(-1, 1, -4),
      xtal::UnitCell(-1, 1, -3),
      xtal::UnitCell(-1, 1, -2),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(-1, 1, 2),
      xtal::UnitCell(-1, 1, 3),
      xtal::UnitCell(-1, 1, 4),
      xtal::UnitCell(-1, 2, -4),
      xtal::UnitCell(-1, 2, -3),
      xtal::UnitCell(-1, 2, -2),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(-1, 2, 1),
      xtal::UnitCell(-1, 2, 2),
      xtal::UnitCell(-1, 2, 3),
      xtal::UnitCell(-1, 3, -4),
      xtal::UnitCell(-1, 3, -3),
      xtal::UnitCell(-1, 3, -2),
      xtal::UnitCell(-1, 3, -1),
      xtal::UnitCell(-1, 3, 0),
      xtal::UnitCell(-1, 3, 1),
      xtal::UnitCell(-1, 3, 2),
      xtal::UnitCell(-1, 4, -4),
      xtal::UnitCell(-1, 4, -3),
      xtal::UnitCell(-1, 4, -2),
      xtal::UnitCell(-1, 4, -1),
      xtal::UnitCell(-1, 4, 0),
      xtal::UnitCell(-1, 4, 1),
      xtal::UnitCell(-1, 5, -3),
      xtal::UnitCell(-1, 5, -2),
      xtal::UnitCell(-1, 5, -1),
      xtal::UnitCell(0, -4, 0),
      xtal::UnitCell(0, -4, 1),
      xtal::UnitCell(0, -4, 2),
      xtal::UnitCell(0, -4, 3),
      xtal::UnitCell(0, -4, 4),
      xtal::UnitCell(0, -3, -1),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, -3, 1),
      xtal::UnitCell(0, -3, 2),
      xtal::UnitCell(0, -3, 3),
      xtal::UnitCell(0, -3, 4),
      xtal::UnitCell(0, -2, -2),
      xtal::UnitCell(0, -2, -1),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, -2, 3),
      xtal::UnitCell(0, -2, 4),
      xtal::UnitCell(0, -1, -3),
      xtal::UnitCell(0, -1, -2),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, -1, 3),
      xtal::UnitCell(0, -1, 4),
      xtal::UnitCell(0, 0, -4),
      xtal::UnitCell(0, 0, -3),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 0, 3),
      xtal::UnitCell(0, 0, 4),
      xtal::UnitCell(0, 1, -4),
      xtal::UnitCell(0, 1, -3),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 1, 2),
      xtal::UnitCell(0, 1, 3),
      xtal::UnitCell(0, 2, -4),
      xtal::UnitCell(0, 2, -3),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(0, 2, 1),
      xtal::UnitCell(0, 2, 2),
      xtal::UnitCell(0, 3, -4),
      xtal::UnitCell(0, 3, -3),
      xtal::UnitCell(0, 3, -2),
      xtal::UnitCell(0, 3, -1),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(0, 3, 1),
      xtal::UnitCell(0, 4, -4),
      xtal::UnitCell(0, 4, -3),
      xtal::UnitCell(0, 4, -2),
      xtal::UnitCell(0, 4, -1),
      xtal::UnitCell(0, 4, 0),
      xtal::UnitCell(1, -5, 1),
      xtal::UnitCell(1, -5, 2),
      xtal::UnitCell(1, -5, 3),
      xtal::UnitCell(1, -4, -1),
      xtal::UnitCell(1, -4, 0),
      xtal::UnitCell(1, -4, 1),
      xtal::UnitCell(1, -4, 2),
      xtal::UnitCell(1, -4, 3),
      xtal::UnitCell(1, -4, 4),
      xtal::UnitCell(1, -3, -2),
      xtal::UnitCell(1, -3, -1),
      xtal::UnitCell(1, -3, 0),
      xtal::UnitCell(1, -3, 1),
      xtal::UnitCell(1, -3, 2),
      xtal::UnitCell(1, -3, 3),
      xtal::UnitCell(1, -3, 4),
      xtal::UnitCell(1, -2, -3),
      xtal::UnitCell(1, -2, -2),
      xtal::UnitCell(1, -2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, -2, 2),
      xtal::UnitCell(1, -2, 3),
      xtal::UnitCell(1, -2, 4),
      xtal::UnitCell(1, -1, -4),
      xtal::UnitCell(1, -1, -3),
      xtal::UnitCell(1, -1, -2),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, -1, 2),
      xtal::UnitCell(1, -1, 3),
      xtal::UnitCell(1, -1, 4),
      xtal::UnitCell(1, 0, -4),
      xtal::UnitCell(1, 0, -3),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 0, 2),
      xtal::UnitCell(1, 0, 3),
      xtal::UnitCell(1, 1, -5),
      xtal::UnitCell(1, 1, -4),
      xtal::UnitCell(1, 1, -3),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, -1),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 1, 1),
      xtal::UnitCell(1, 1, 2),
      xtal::UnitCell(1, 1, 3),
      xtal::UnitCell(1, 2, -5),
      xtal::UnitCell(1, 2, -4),
      xtal::UnitCell(1, 2, -3),
      xtal::UnitCell(1, 2, -2),
      xtal::UnitCell(1, 2, -1),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(1, 2, 1),
      xtal::UnitCell(1, 2, 2),
      xtal::UnitCell(1, 3, -5),
      xtal::UnitCell(1, 3, -4),
      xtal::UnitCell(1, 3, -3),
      xtal::UnitCell(1, 3, -2),
      xtal::UnitCell(1, 3, -1),
      xtal::UnitCell(1, 3, 0),
      xtal::UnitCell(1, 3, 1),
      xtal::UnitCell(1, 4, -4),
      xtal::UnitCell(1, 4, -3),
      xtal::UnitCell(1, 4, -2),
      xtal::UnitCell(1, 4, -1),
      xtal::UnitCell(2, -5, 1),
      xtal::UnitCell(2, -5, 2),
      xtal::UnitCell(2, -4, -1),
      xtal::UnitCell(2, -4, 0),
      xtal::UnitCell(2, -4, 1),
      xtal::UnitCell(2, -4, 2),
      xtal::UnitCell(2, -4, 3),
      xtal::UnitCell(2, -3, -2),
      xtal::UnitCell(2, -3, -1),
      xtal::UnitCell(2, -3, 0),
      xtal::UnitCell(2, -3, 1),
      xtal::UnitCell(2, -3, 2),
      xtal::UnitCell(2, -3, 3),
      xtal::UnitCell(2, -2, -3),
      xtal::UnitCell(2, -2, -2),
      xtal::UnitCell(2, -2, -1),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, -2, 1),
      xtal::UnitCell(2, -2, 2),
      xtal::UnitCell(2, -2, 3),
      xtal::UnitCell(2, -1, -4),
      xtal::UnitCell(2, -1, -3),
      xtal::UnitCell(2, -1, -2),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, -1, 1),
      xtal::UnitCell(2, -1, 2),
      xtal::UnitCell(2, -1, 3),
      xtal::UnitCell(2, 0, -4),
      xtal::UnitCell(2, 0, -3),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, -1),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 0, 1),
      xtal::UnitCell(2, 0, 2),
      xtal::UnitCell(2, 1, -5),
      xtal::UnitCell(2, 1, -4),
      xtal::UnitCell(2, 1, -3),
      xtal::UnitCell(2, 1, -2),
      xtal::UnitCell(2, 1, -1),
      xtal::UnitCell(2, 1, 0),
      xtal::UnitCell(2, 1, 1),
      xtal::UnitCell(2, 1, 2),
      xtal::UnitCell(2, 2, -5),
      xtal::UnitCell(2, 2, -4),
      xtal::UnitCell(2, 2, -3),
      xtal::UnitCell(2, 2, -2),
      xtal::UnitCell(2, 2, -1),
      xtal::UnitCell(2, 2, 0),
      xtal::UnitCell(2, 2, 1),
      xtal::UnitCell(2, 3, -4),
      xtal::UnitCell(2, 3, -3),
      xtal::UnitCell(2, 3, -2),
      xtal::UnitCell(2, 3, -1),
      xtal::UnitCell(3, -5, 1),
      xtal::UnitCell(3, -4, -1),
      xtal::UnitCell(3, -4, 0),
      xtal::UnitCell(3, -4, 1),
      xtal::UnitCell(3, -4, 2),
      xtal::UnitCell(3, -3, -3),
      xtal::UnitCell(3, -3, -2),
      xtal::UnitCell(3, -3, -1),
      xtal::UnitCell(3, -3, 0),
      xtal::UnitCell(3, -3, 1),
      xtal::UnitCell(3, -3, 2),
      xtal::UnitCell(3, -3, 3),
      xtal::UnitCell(3, -2, -3),
      xtal::UnitCell(3, -2, -2),
      xtal::UnitCell(3, -2, -1),
      xtal::UnitCell(3, -2, 0),
      xtal::UnitCell(3, -2, 1),
      xtal::UnitCell(3, -2, 2),
      xtal::UnitCell(3, -1, -4),
      xtal::UnitCell(3, -1, -3),
      xtal::UnitCell(3, -1, -2),
      xtal::UnitCell(3, -1, -1),
      xtal::UnitCell(3, -1, 0),
      xtal::UnitCell(3, -1, 1),
      xtal::UnitCell(3, -1, 2),
      xtal::UnitCell(3, 0, -4),
      xtal::UnitCell(3, 0, -3),
      xtal::UnitCell(3, 0, -2),
      xtal::UnitCell(3, 0, -1),
      xtal::UnitCell(3, 0, 0),
      xtal::UnitCell(3, 0, 1),
      xtal::UnitCell(3, 1, -5),
      xtal::UnitCell(3, 1, -4),
      xtal::UnitCell(3, 1, -3),
      xtal::UnitCell(3, 1, -2),
      xtal::UnitCell(3, 1, -1),
      xtal::UnitCell(3, 1, 0),
      xtal::UnitCell(3, 1, 1),
      xtal::UnitCell(3, 2, -4),
      xtal::UnitCell(3, 2, -3),
      xtal::UnitCell(3, 2, -2),
      xtal::UnitCell(3, 2, -1),
      xtal::UnitCell(3, 3, -3),
      xtal::UnitCell(4, -4, -1),
      xtal::UnitCell(4, -4, 0),
      xtal::UnitCell(4, -4, 1),
      xtal::UnitCell(4, -3, -2),
      xtal::UnitCell(4, -3, -1),
      xtal::UnitCell(4, -3, 0),
      xtal::UnitCell(4, -3, 1),
      xtal::UnitCell(4, -2, -3),
      xtal::UnitCell(4, -2, -2),
      xtal::UnitCell(4, -2, -1),
      xtal::UnitCell(4, -2, 0),
      xtal::UnitCell(4, -2, 1),
      xtal::UnitCell(4, -1, -4),
      xtal::UnitCell(4, -1, -3),
      xtal::UnitCell(4, -1, -2),
      xtal::UnitCell(4, -1, -1),
      xtal::UnitCell(4, -1, 0),
      xtal::UnitCell(4, -1, 1),
      xtal::UnitCell(4, 0, -4),
      xtal::UnitCell(4, 0, -3),
      xtal::UnitCell(4, 0, -2),
      xtal::UnitCell(4, 0, -1),
      xtal::UnitCell(4, 0, 0),
      xtal::UnitCell(4, 1, -4),
      xtal::UnitCell(4, 1, -3),
      xtal::UnitCell(4, 1, -2),
      xtal::UnitCell(4, 1, -1),
      xtal::UnitCell(5, -3, -1),
      xtal::UnitCell(5, -2, -2),
      xtal::UnitCell(5, -2, -1),
      xtal::UnitCell(5, -1, -3),
      xtal::UnitCell(5, -1, -2),
      xtal::UnitCell(5, -1, -1)
    };


    m_orbit_neighborhood.resize(corr_size());
    m_orbit_site_neighborhood.resize(corr_size());
    m_orbit_neighborhood[1] = std::set<xtal::UnitCell> {
      xtal::UnitCell(0, 0, 0)
    };

    m_orbit_site_neighborhood[1] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, 0, 0, 0)
    };

    m_orbit_neighborhood[2] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0)
    };

    m_orbit_site_neighborhood[2] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0)
    };

    m_orbit_neighborhood[3] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 1, -1)
    };

    m_orbit_site_neighborhood[3] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 1, -1)
    };

    m_orbit_neighborhood[4] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -1)
    };

    m_orbit_site_neighborhood[4] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -1)
    };

    m_orbit_neighborhood[5] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[5] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 2),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -2, 2),
      xtal::UnitCellCoord(0, 0, 0, -2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 2),
      xtal::UnitCellCoord(0, 0, 2, -2),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 0, -2),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[6] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 1),
      xtal::UnitCell(-2, -1, 2),
      xtal::UnitCell(-2, 1, -1),
      xtal::UnitCell(-2, 1, 2),
      xtal::UnitCell(-2, 2, -1),
      xtal::UnitCell(-2, 2, 1),
      xtal::UnitCell(-1, -2, 1),
      xtal::UnitCell(-1, -2, 2),
      xtal::UnitCell(-1, 1, -2),
      xtal::UnitCell(-1, 1, 2),
      xtal::UnitCell(-1, 2, -2),
      xtal::UnitCell(-1, 2, 1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -2, -1),
      xtal::UnitCell(1, -2, 2),
      xtal::UnitCell(1, -1, -2),
      xtal::UnitCell(1, -1, 2),
      xtal::UnitCell(1, 2, -2),
      xtal::UnitCell(1, 2, -1),
      xtal::UnitCell(2, -2, -1),
      xtal::UnitCell(2, -2, 1),
      xtal::UnitCell(2, -1, -2),
      xtal::UnitCell(2, -1, 1),
      xtal::UnitCell(2, 1, -2),
      xtal::UnitCell(2, 1, -1)
    };

    m_orbit_site_neighborhood[6] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 1),
      xtal::UnitCellCoord(0, -2, -1, 2),
      xtal::UnitCellCoord(0, -2, 1, -1),
      xtal::UnitCellCoord(0, -2, 1, 2),
      xtal::UnitCellCoord(0, -2, 2, -1),
      xtal::UnitCellCoord(0, -2, 2, 1),
      xtal::UnitCellCoord(0, -1, -2, 1),
      xtal::UnitCellCoord(0, -1, -2, 2),
      xtal::UnitCellCoord(0, -1, 1, -2),
      xtal::UnitCellCoord(0, -1, 1, 2),
      xtal::UnitCellCoord(0, -1, 2, -2),
      xtal::UnitCellCoord(0, -1, 2, 1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -2, -1),
      xtal::UnitCellCoord(0, 1, -2, 2),
      xtal::UnitCellCoord(0, 1, -1, -2),
      xtal::UnitCellCoord(0, 1, -1, 2),
      xtal::UnitCellCoord(0, 1, 2, -2),
      xtal::UnitCellCoord(0, 1, 2, -1),
      xtal::UnitCellCoord(0, 2, -2, -1),
      xtal::UnitCellCoord(0, 2, -2, 1),
      xtal::UnitCellCoord(0, 2, -1, -2),
      xtal::UnitCellCoord(0, 2, -1, 1),
      xtal::UnitCellCoord(0, 2, 1, -2),
      xtal::UnitCellCoord(0, 2, 1, -1)
    };

    m_orbit_neighborhood[7] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 1, 1),
      xtal::UnitCell(-1, -1, -1),
      xtal::UnitCell(-1, -1, 3),
      xtal::UnitCell(-1, 3, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -3, 1),
      xtal::UnitCell(1, 1, -3),
      xtal::UnitCell(1, 1, 1),
      xtal::UnitCell(3, -1, -1)
    };

    m_orbit_site_neighborhood[7] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, -1),
      xtal::UnitCellCoord(0, -1, -1, 3),
      xtal::UnitCellCoord(0, -1, 3, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -3, 1),
      xtal::UnitCellCoord(0, 1, 1, -3),
      xtal::UnitCellCoord(0, 1, 1, 1),
      xtal::UnitCellCoord(0, 3, -1, -1)
    };

    m_orbit_neighborhood[8] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 1),
      xtal::UnitCell(-3, 0, 2),
      xtal::UnitCell(-3, 1, 0),
      xtal::UnitCell(-3, 1, 2),
      xtal::UnitCell(-3, 2, 0),
      xtal::UnitCell(-3, 2, 1),
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, -1, 3),
      xtal::UnitCell(-2, 0, -1),
      xtal::UnitCell(-2, 0, 3),
      xtal::UnitCell(-2, 3, -1),
      xtal::UnitCell(-2, 3, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -2, 3),
      xtal::UnitCell(-1, 0, -2),
      xtal::UnitCell(-1, 0, 3),
      xtal::UnitCell(-1, 3, -2),
      xtal::UnitCell(-1, 3, 0),
      xtal::UnitCell(0, -3, 1),
      xtal::UnitCell(0, -3, 2),
      xtal::UnitCell(0, -2, -1),
      xtal::UnitCell(0, -2, 3),
      xtal::UnitCell(0, -1, -2),
      xtal::UnitCell(0, -1, 3),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, -3),
      xtal::UnitCell(0, 1, 2),
      xtal::UnitCell(0, 2, -3),
      xtal::UnitCell(0, 2, 1),
      xtal::UnitCell(0, 3, -2),
      xtal::UnitCell(0, 3, -1),
      xtal::UnitCell(1, -3, 0),
      xtal::UnitCell(1, -3, 2),
      xtal::UnitCell(1, 0, -3),
      xtal::UnitCell(1, 0, 2),
      xtal::UnitCell(1, 2, -3),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -3, 0),
      xtal::UnitCell(2, -3, 1),
      xtal::UnitCell(2, 0, -3),
      xtal::UnitCell(2, 0, 1),
      xtal::UnitCell(2, 1, -3),
      xtal::UnitCell(2, 1, 0),
      xtal::UnitCell(3, -2, -1),
      xtal::UnitCell(3, -2, 0),
      xtal::UnitCell(3, -1, -2),
      xtal::UnitCell(3, -1, 0),
      xtal::UnitCell(3, 0, -2),
      xtal::UnitCell(3, 0, -1)
    };

    m_orbit_site_neighborhood[8] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 0, 1),
      xtal::UnitCellCoord(0, -3, 0, 2),
      xtal::UnitCellCoord(0, -3, 1, 0),
      xtal::UnitCellCoord(0, -3, 1, 2),
      xtal::UnitCellCoord(0, -3, 2, 0),
      xtal::UnitCellCoord(0, -3, 2, 1),
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, -1, 3),
      xtal::UnitCellCoord(0, -2, 0, -1),
      xtal::UnitCellCoord(0, -2, 0, 3),
      xtal::UnitCellCoord(0, -2, 3, -1),
      xtal::UnitCellCoord(0, -2, 3, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -2, 3),
      xtal::UnitCellCoord(0, -1, 0, -2),
      xtal::UnitCellCoord(0, -1, 0, 3),
      xtal::UnitCellCoord(0, -1, 3, -2),
      xtal::UnitCellCoord(0, -1, 3, 0),
      xtal::UnitCellCoord(0, 0, -3, 1),
      xtal::UnitCellCoord(0, 0, -3, 2),
      xtal::UnitCellCoord(0, 0, -2, -1),
      xtal::UnitCellCoord(0, 0, -2, 3),
      xtal::UnitCellCoord(0, 0, -1, -2),
      xtal::UnitCellCoord(0, 0, -1, 3),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, -3),
      xtal::UnitCellCoord(0, 0, 1, 2),
      xtal::UnitCellCoord(0, 0, 2, -3),
      xtal::UnitCellCoord(0, 0, 2, 1),
      xtal::UnitCellCoord(0, 0, 3, -2),
      xtal::UnitCellCoord(0, 0, 3, -1),
      xtal::UnitCellCoord(0, 1, -3, 0),
      xtal::UnitCellCoord(0, 1, -3, 2),
      xtal::UnitCellCoord(0, 1, 0, -3),
      xtal::UnitCellCoord(0, 1, 0, 2),
      xtal::UnitCellCoord(0, 1, 2, -3),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -3, 0),
      xtal::UnitCellCoord(0, 2, -3, 1),
      xtal::UnitCellCoord(0, 2, 0, -3),
      xtal::UnitCellCoord(0, 2, 0, 1),
      xtal::UnitCellCoord(0, 2, 1, -3),
      xtal::UnitCellCoord(0, 2, 1, 0),
      xtal::UnitCellCoord(0, 3, -2, -1),
      xtal::UnitCellCoord(0, 3, -2, 0),
      xtal::UnitCellCoord(0, 3, -1, -2),
      xtal::UnitCellCoord(0, 3, -1, 0),
      xtal::UnitCellCoord(0, 3, 0, -2),
      xtal::UnitCellCoord(0, 3, 0, -1)
    };

    m_orbit_neighborhood[9] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -2, 2),
      xtal::UnitCell(-2, 2, -2),
      xtal::UnitCell(-2, 2, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(2, -2, -2),
      xtal::UnitCell(2, -2, 2),
      xtal::UnitCell(2, 2, -2)
    };

    m_orbit_site_neighborhood[9] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -2, 2),
      xtal::UnitCellCoord(0, -2, 2, -2),
      xtal::UnitCellCoord(0, -2, 2, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 2, -2, -2),
      xtal::UnitCellCoord(0, 2, -2, 2),
      xtal::UnitCellCoord(0, 2, 2, -2)
    };

    m_orbit_neighborhood[10] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(-3, 0, 3),
      xtal::UnitCell(-3, 3, 0),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, -3, 3),
      xtal::UnitCell(0, 0, -3),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 3),
      xtal::UnitCell(0, 3, -3),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(3, -3, 0),
      xtal::UnitCell(3, 0, -3),
      xtal::UnitCell(3, 0, 0)
    };

    m_orbit_site_neighborhood[10] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 0, 0),
      xtal::UnitCellCoord(0, -3, 0, 3),
      xtal::UnitCellCoord(0, -3, 3, 0),
      xtal::UnitCellCoord(0, 0, -3, 0),
      xtal::UnitCellCoord(0, 0, -3, 3),
      xtal::UnitCellCoord(0, 0, 0, -3),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 3),
      xtal::UnitCellCoord(0, 0, 3, -3),
      xtal::UnitCellCoord(0, 0, 3, 0),
      xtal::UnitCellCoord(0, 3, -3, 0),
      xtal::UnitCellCoord(0, 3, 0, -3),
      xtal::UnitCellCoord(0, 3, 0, 0)
    };

    m_orbit_neighborhood[11] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, -1, 2),
      xtal::UnitCell(-3, 2, -1),
      xtal::UnitCell(-3, 2, 2),
      xtal::UnitCell(-2, -2, 1),
      xtal::UnitCell(-2, -2, 3),
      xtal::UnitCell(-2, 1, -2),
      xtal::UnitCell(-2, 1, 3),
      xtal::UnitCell(-2, 3, -2),
      xtal::UnitCell(-2, 3, 1),
      xtal::UnitCell(-1, -3, 2),
      xtal::UnitCell(-1, 2, -3),
      xtal::UnitCell(-1, 2, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -2, -2),
      xtal::UnitCell(1, -2, 3),
      xtal::UnitCell(1, 3, -2),
      xtal::UnitCell(2, -3, -1),
      xtal::UnitCell(2, -3, 2),
      xtal::UnitCell(2, -1, -3),
      xtal::UnitCell(2, -1, 2),
      xtal::UnitCell(2, 2, -3),
      xtal::UnitCell(2, 2, -1),
      xtal::UnitCell(3, -2, -2),
      xtal::UnitCell(3, -2, 1),
      xtal::UnitCell(3, 1, -2)
    };

    m_orbit_site_neighborhood[11] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, -1, 2),
      xtal::UnitCellCoord(0, -3, 2, -1),
      xtal::UnitCellCoord(0, -3, 2, 2),
      xtal::UnitCellCoord(0, -2, -2, 1),
      xtal::UnitCellCoord(0, -2, -2, 3),
      xtal::UnitCellCoord(0, -2, 1, -2),
      xtal::UnitCellCoord(0, -2, 1, 3),
      xtal::UnitCellCoord(0, -2, 3, -2),
      xtal::UnitCellCoord(0, -2, 3, 1),
      xtal::UnitCellCoord(0, -1, -3, 2),
      xtal::UnitCellCoord(0, -1, 2, -3),
      xtal::UnitCellCoord(0, -1, 2, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -2, -2),
      xtal::UnitCellCoord(0, 1, -2, 3),
      xtal::UnitCellCoord(0, 1, 3, -2),
      xtal::UnitCellCoord(0, 2, -3, -1),
      xtal::UnitCellCoord(0, 2, -3, 2),
      xtal::UnitCellCoord(0, 2, -1, -3),
      xtal::UnitCellCoord(0, 2, -1, 2),
      xtal::UnitCellCoord(0, 2, 2, -3),
      xtal::UnitCellCoord(0, 2, 2, -1),
      xtal::UnitCellCoord(0, 3, -2, -2),
      xtal::UnitCellCoord(0, 3, -2, 1),
      xtal::UnitCellCoord(0, 3, 1, -2)
    };

    m_orbit_neighborhood[12] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, -1, 1),
      xtal::UnitCell(-3, -1, 3),
      xtal::UnitCell(-3, 1, -1),
      xtal::UnitCell(-3, 1, 3),
      xtal::UnitCell(-3, 3, -1),
      xtal::UnitCell(-3, 3, 1),
      xtal::UnitCell(-1, -3, 1),
      xtal::UnitCell(-1, -3, 3),
      xtal::UnitCell(-1, 1, -3),
      xtal::UnitCell(-1, 1, 3),
      xtal::UnitCell(-1, 3, -3),
      xtal::UnitCell(-1, 3, 1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -3, -1),
      xtal::UnitCell(1, -3, 3),
      xtal::UnitCell(1, -1, -3),
      xtal::UnitCell(1, -1, 3),
      xtal::UnitCell(1, 3, -3),
      xtal::UnitCell(1, 3, -1),
      xtal::UnitCell(3, -3, -1),
      xtal::UnitCell(3, -3, 1),
      xtal::UnitCell(3, -1, -3),
      xtal::UnitCell(3, -1, 1),
      xtal::UnitCell(3, 1, -3),
      xtal::UnitCell(3, 1, -1)
    };

    m_orbit_site_neighborhood[12] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, -1, 1),
      xtal::UnitCellCoord(0, -3, -1, 3),
      xtal::UnitCellCoord(0, -3, 1, -1),
      xtal::UnitCellCoord(0, -3, 1, 3),
      xtal::UnitCellCoord(0, -3, 3, -1),
      xtal::UnitCellCoord(0, -3, 3, 1),
      xtal::UnitCellCoord(0, -1, -3, 1),
      xtal::UnitCellCoord(0, -1, -3, 3),
      xtal::UnitCellCoord(0, -1, 1, -3),
      xtal::UnitCellCoord(0, -1, 1, 3),
      xtal::UnitCellCoord(0, -1, 3, -3),
      xtal::UnitCellCoord(0, -1, 3, 1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -3, -1),
      xtal::UnitCellCoord(0, 1, -3, 3),
      xtal::UnitCellCoord(0, 1, -1, -3),
      xtal::UnitCellCoord(0, 1, -1, 3),
      xtal::UnitCellCoord(0, 1, 3, -3),
      xtal::UnitCellCoord(0, 1, 3, -1),
      xtal::UnitCellCoord(0, 3, -3, -1),
      xtal::UnitCellCoord(0, 3, -3, 1),
      xtal::UnitCellCoord(0, 3, -1, -3),
      xtal::UnitCellCoord(0, 3, -1, 1),
      xtal::UnitCellCoord(0, 3, 1, -3),
      xtal::UnitCellCoord(0, 3, 1, -1)
    };

    m_orbit_neighborhood[13] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-4, 1, 1),
      xtal::UnitCell(-4, 1, 2),
      xtal::UnitCell(-4, 2, 1),
      xtal::UnitCell(-2, -1, -1),
      xtal::UnitCell(-2, -1, 4),
      xtal::UnitCell(-2, 4, -1),
      xtal::UnitCell(-1, -2, -1),
      xtal::UnitCell(-1, -2, 4),
      xtal::UnitCell(-1, -1, -2),
      xtal::UnitCell(-1, -1, 4),
      xtal::UnitCell(-1, 4, -2),
      xtal::UnitCell(-1, 4, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -4, 1),
      xtal::UnitCell(1, -4, 2),
      xtal::UnitCell(1, 1, -4),
      xtal::UnitCell(1, 1, 2),
      xtal::UnitCell(1, 2, -4),
      xtal::UnitCell(1, 2, 1),
      xtal::UnitCell(2, -4, 1),
      xtal::UnitCell(2, 1, -4),
      xtal::UnitCell(2, 1, 1),
      xtal::UnitCell(4, -2, -1),
      xtal::UnitCell(4, -1, -2),
      xtal::UnitCell(4, -1, -1)
    };

    m_orbit_site_neighborhood[13] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -4, 1, 1),
      xtal::UnitCellCoord(0, -4, 1, 2),
      xtal::UnitCellCoord(0, -4, 2, 1),
      xtal::UnitCellCoord(0, -2, -1, -1),
      xtal::UnitCellCoord(0, -2, -1, 4),
      xtal::UnitCellCoord(0, -2, 4, -1),
      xtal::UnitCellCoord(0, -1, -2, -1),
      xtal::UnitCellCoord(0, -1, -2, 4),
      xtal::UnitCellCoord(0, -1, -1, -2),
      xtal::UnitCellCoord(0, -1, -1, 4),
      xtal::UnitCellCoord(0, -1, 4, -2),
      xtal::UnitCellCoord(0, -1, 4, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -4, 1),
      xtal::UnitCellCoord(0, 1, -4, 2),
      xtal::UnitCellCoord(0, 1, 1, -4),
      xtal::UnitCellCoord(0, 1, 1, 2),
      xtal::UnitCellCoord(0, 1, 2, -4),
      xtal::UnitCellCoord(0, 1, 2, 1),
      xtal::UnitCellCoord(0, 2, -4, 1),
      xtal::UnitCellCoord(0, 2, 1, -4),
      xtal::UnitCellCoord(0, 2, 1, 1),
      xtal::UnitCellCoord(0, 4, -2, -1),
      xtal::UnitCellCoord(0, 4, -1, -2),
      xtal::UnitCellCoord(0, 4, -1, -1)
    };

    m_orbit_neighborhood[14] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-4, 0, 2),
      xtal::UnitCell(-4, 2, 0),
      xtal::UnitCell(-4, 2, 2),
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, -2, 4),
      xtal::UnitCell(-2, 0, -2),
      xtal::UnitCell(-2, 0, 4),
      xtal::UnitCell(-2, 4, -2),
      xtal::UnitCell(-2, 4, 0),
      xtal::UnitCell(0, -4, 2),
      xtal::UnitCell(0, -2, -2),
      xtal::UnitCell(0, -2, 4),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 2, -4),
      xtal::UnitCell(0, 2, 2),
      xtal::UnitCell(0, 4, -2),
      xtal::UnitCell(2, -4, 0),
      xtal::UnitCell(2, -4, 2),
      xtal::UnitCell(2, 0, -4),
      xtal::UnitCell(2, 0, 2),
      xtal::UnitCell(2, 2, -4),
      xtal::UnitCell(2, 2, 0),
      xtal::UnitCell(4, -2, -2),
      xtal::UnitCell(4, -2, 0),
      xtal::UnitCell(4, 0, -2)
    };

    m_orbit_site_neighborhood[14] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -4, 0, 2),
      xtal::UnitCellCoord(0, -4, 2, 0),
      xtal::UnitCellCoord(0, -4, 2, 2),
      xtal::UnitCellCoord(0, -2, -2, 0),
      xtal::UnitCellCoord(0, -2, -2, 4),
      xtal::UnitCellCoord(0, -2, 0, -2),
      xtal::UnitCellCoord(0, -2, 0, 4),
      xtal::UnitCellCoord(0, -2, 4, -2),
      xtal::UnitCellCoord(0, -2, 4, 0),
      xtal::UnitCellCoord(0, 0, -4, 2),
      xtal::UnitCellCoord(0, 0, -2, -2),
      xtal::UnitCellCoord(0, 0, -2, 4),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 2, -4),
      xtal::UnitCellCoord(0, 0, 2, 2),
      xtal::UnitCellCoord(0, 0, 4, -2),
      xtal::UnitCellCoord(0, 2, -4, 0),
      xtal::UnitCellCoord(0, 2, -4, 2),
      xtal::UnitCellCoord(0, 2, 0, -4),
      xtal::UnitCellCoord(0, 2, 0, 2),
      xtal::UnitCellCoord(0, 2, 2, -4),
      xtal::UnitCellCoord(0, 2, 2, 0),
      xtal::UnitCellCoord(0, 4, -2, -2),
      xtal::UnitCellCoord(0, 4, -2, 0),
      xtal::UnitCellCoord(0, 4, 0, -2)
    };

    m_orbit_neighborhood[15] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, -2, 2),
      xtal::UnitCell(-3, -2, 3),
      xtal::UnitCell(-3, 2, -2),
      xtal::UnitCell(-3, 2, 3),
      xtal::UnitCell(-3, 3, -2),
      xtal::UnitCell(-3, 3, 2),
      xtal::UnitCell(-2, -3, 2),
      xtal::UnitCell(-2, -3, 3),
      xtal::UnitCell(-2, 2, -3),
      xtal::UnitCell(-2, 2, 3),
      xtal::UnitCell(-2, 3, -3),
      xtal::UnitCell(-2, 3, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(2, -3, -2),
      xtal::UnitCell(2, -3, 3),
      xtal::UnitCell(2, -2, -3),
      xtal::UnitCell(2, -2, 3),
      xtal::UnitCell(2, 3, -3),
      xtal::UnitCell(2, 3, -2),
      xtal::UnitCell(3, -3, -2),
      xtal::UnitCell(3, -3, 2),
      xtal::UnitCell(3, -2, -3),
      xtal::UnitCell(3, -2, 2),
      xtal::UnitCell(3, 2, -3),
      xtal::UnitCell(3, 2, -2)
    };

    m_orbit_site_neighborhood[15] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, -2, 2),
      xtal::UnitCellCoord(0, -3, -2, 3),
      xtal::UnitCellCoord(0, -3, 2, -2),
      xtal::UnitCellCoord(0, -3, 2, 3),
      xtal::UnitCellCoord(0, -3, 3, -2),
      xtal::UnitCellCoord(0, -3, 3, 2),
      xtal::UnitCellCoord(0, -2, -3, 2),
      xtal::UnitCellCoord(0, -2, -3, 3),
      xtal::UnitCellCoord(0, -2, 2, -3),
      xtal::UnitCellCoord(0, -2, 2, 3),
      xtal::UnitCellCoord(0, -2, 3, -3),
      xtal::UnitCellCoord(0, -2, 3, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 2, -3, -2),
      xtal::UnitCellCoord(0, 2, -3, 3),
      xtal::UnitCellCoord(0, 2, -2, -3),
      xtal::UnitCellCoord(0, 2, -2, 3),
      xtal::UnitCellCoord(0, 2, 3, -3),
      xtal::UnitCellCoord(0, 2, 3, -2),
      xtal::UnitCellCoord(0, 3, -3, -2),
      xtal::UnitCellCoord(0, 3, -3, 2),
      xtal::UnitCellCoord(0, 3, -2, -3),
      xtal::UnitCellCoord(0, 3, -2, 2),
      xtal::UnitCellCoord(0, 3, 2, -3),
      xtal::UnitCellCoord(0, 3, 2, -2)
    };

    m_orbit_neighborhood[16] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-4, 0, 1),
      xtal::UnitCell(-4, 0, 3),
      xtal::UnitCell(-4, 1, 0),
      xtal::UnitCell(-4, 1, 3),
      xtal::UnitCell(-4, 3, 0),
      xtal::UnitCell(-4, 3, 1),
      xtal::UnitCell(-3, -1, 0),
      xtal::UnitCell(-3, -1, 4),
      xtal::UnitCell(-3, 0, -1),
      xtal::UnitCell(-3, 0, 4),
      xtal::UnitCell(-3, 4, -1),
      xtal::UnitCell(-3, 4, 0),
      xtal::UnitCell(-1, -3, 0),
      xtal::UnitCell(-1, -3, 4),
      xtal::UnitCell(-1, 0, -3),
      xtal::UnitCell(-1, 0, 4),
      xtal::UnitCell(-1, 4, -3),
      xtal::UnitCell(-1, 4, 0),
      xtal::UnitCell(0, -4, 1),
      xtal::UnitCell(0, -4, 3),
      xtal::UnitCell(0, -3, -1),
      xtal::UnitCell(0, -3, 4),
      xtal::UnitCell(0, -1, -3),
      xtal::UnitCell(0, -1, 4),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, -4),
      xtal::UnitCell(0, 1, 3),
      xtal::UnitCell(0, 3, -4),
      xtal::UnitCell(0, 3, 1),
      xtal::UnitCell(0, 4, -3),
      xtal::UnitCell(0, 4, -1),
      xtal::UnitCell(1, -4, 0),
      xtal::UnitCell(1, -4, 3),
      xtal::UnitCell(1, 0, -4),
      xtal::UnitCell(1, 0, 3),
      xtal::UnitCell(1, 3, -4),
      xtal::UnitCell(1, 3, 0),
      xtal::UnitCell(3, -4, 0),
      xtal::UnitCell(3, -4, 1),
      xtal::UnitCell(3, 0, -4),
      xtal::UnitCell(3, 0, 1),
      xtal::UnitCell(3, 1, -4),
      xtal::UnitCell(3, 1, 0),
      xtal::UnitCell(4, -3, -1),
      xtal::UnitCell(4, -3, 0),
      xtal::UnitCell(4, -1, -3),
      xtal::UnitCell(4, -1, 0),
      xtal::UnitCell(4, 0, -3),
      xtal::UnitCell(4, 0, -1)
    };

    m_orbit_site_neighborhood[16] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -4, 0, 1),
      xtal::UnitCellCoord(0, -4, 0, 3),
      xtal::UnitCellCoord(0, -4, 1, 0),
      xtal::UnitCellCoord(0, -4, 1, 3),
      xtal::UnitCellCoord(0, -4, 3, 0),
      xtal::UnitCellCoord(0, -4, 3, 1),
      xtal::UnitCellCoord(0, -3, -1, 0),
      xtal::UnitCellCoord(0, -3, -1, 4),
      xtal::UnitCellCoord(0, -3, 0, -1),
      xtal::UnitCellCoord(0, -3, 0, 4),
      xtal::UnitCellCoord(0, -3, 4, -1),
      xtal::UnitCellCoord(0, -3, 4, 0),
      xtal::UnitCellCoord(0, -1, -3, 0),
      xtal::UnitCellCoord(0, -1, -3, 4),
      xtal::UnitCellCoord(0, -1, 0, -3),
      xtal::UnitCellCoord(0, -1, 0, 4),
      xtal::UnitCellCoord(0, -1, 4, -3),
      xtal::UnitCellCoord(0, -1, 4, 0),
      xtal::UnitCellCoord(0, 0, -4, 1),
      xtal::UnitCellCoord(0, 0, -4, 3),
      xtal::UnitCellCoord(0, 0, -3, -1),
      xtal::UnitCellCoord(0, 0, -3, 4),
      xtal::UnitCellCoord(0, 0, -1, -3),
      xtal::UnitCellCoord(0, 0, -1, 4),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, -4),
      xtal::UnitCellCoord(0, 0, 1, 3),
      xtal::UnitCellCoord(0, 0, 3, -4),
      xtal::UnitCellCoord(0, 0, 3, 1),
      xtal::UnitCellCoord(0, 0, 4, -3),
      xtal::UnitCellCoord(0, 0, 4, -1),
      xtal::UnitCellCoord(0, 1, -4, 0),
      xtal::UnitCellCoord(0, 1, -4, 3),
      xtal::UnitCellCoord(0, 1, 0, -4),
      xtal::UnitCellCoord(0, 1, 0, 3),
      xtal::UnitCellCoord(0, 1, 3, -4),
      xtal::UnitCellCoord(0, 1, 3, 0),
      xtal::UnitCellCoord(0, 3, -4, 0),
      xtal::UnitCellCoord(0, 3, -4, 1),
      xtal::UnitCellCoord(0, 3, 0, -4),
      xtal::UnitCellCoord(0, 3, 0, 1),
      xtal::UnitCellCoord(0, 3, 1, -4),
      xtal::UnitCellCoord(0, 3, 1, 0),
      xtal::UnitCellCoord(0, 4, -3, -1),
      xtal::UnitCellCoord(0, 4, -3, 0),
      xtal::UnitCellCoord(0, 4, -1, -3),
      xtal::UnitCellCoord(0, 4, -1, 0),
      xtal::UnitCellCoord(0, 4, 0, -3),
      xtal::UnitCellCoord(0, 4, 0, -1)
    };

    m_orbit_neighborhood[17] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-4, -1, 2),
      xtal::UnitCell(-4, -1, 3),
      xtal::UnitCell(-4, 2, -1),
      xtal::UnitCell(-4, 2, 3),
      xtal::UnitCell(-4, 3, -1),
      xtal::UnitCell(-4, 3, 2),
      xtal::UnitCell(-3, -2, 1),
      xtal::UnitCell(-3, -2, 4),
      xtal::UnitCell(-3, 1, -2),
      xtal::UnitCell(-3, 1, 4),
      xtal::UnitCell(-3, 4, -2),
      xtal::UnitCell(-3, 4, 1),
      xtal::UnitCell(-2, -3, 1),
      xtal::UnitCell(-2, -3, 4),
      xtal::UnitCell(-2, 1, -3),
      xtal::UnitCell(-2, 1, 4),
      xtal::UnitCell(-2, 4, -3),
      xtal::UnitCell(-2, 4, 1),
      xtal::UnitCell(-1, -4, 2),
      xtal::UnitCell(-1, -4, 3),
      xtal::UnitCell(-1, 2, -4),
      xtal::UnitCell(-1, 2, 3),
      xtal::UnitCell(-1, 3, -4),
      xtal::UnitCell(-1, 3, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -3, -2),
      xtal::UnitCell(1, -3, 4),
      xtal::UnitCell(1, -2, -3),
      xtal::UnitCell(1, -2, 4),
      xtal::UnitCell(1, 4, -3),
      xtal::UnitCell(1, 4, -2),
      xtal::UnitCell(2, -4, -1),
      xtal::UnitCell(2, -4, 3),
      xtal::UnitCell(2, -1, -4),
      xtal::UnitCell(2, -1, 3),
      xtal::UnitCell(2, 3, -4),
      xtal::UnitCell(2, 3, -1),
      xtal::UnitCell(3, -4, -1),
      xtal::UnitCell(3, -4, 2),
      xtal::UnitCell(3, -1, -4),
      xtal::UnitCell(3, -1, 2),
      xtal::UnitCell(3, 2, -4),
      xtal::UnitCell(3, 2, -1),
      xtal::UnitCell(4, -3, -2),
      xtal::UnitCell(4, -3, 1),
      xtal::UnitCell(4, -2, -3),
      xtal::UnitCell(4, -2, 1),
      xtal::UnitCell(4, 1, -3),
      xtal::UnitCell(4, 1, -2)
    };

    m_orbit_site_neighborhood[17] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -4, -1, 2),
      xtal::UnitCellCoord(0, -4, -1, 3),
      xtal::UnitCellCoord(0, -4, 2, -1),
      xtal::UnitCellCoord(0, -4, 2, 3),
      xtal::UnitCellCoord(0, -4, 3, -1),
      xtal::UnitCellCoord(0, -4, 3, 2),
      xtal::UnitCellCoord(0, -3, -2, 1),
      xtal::UnitCellCoord(0, -3, -2, 4),
      xtal::UnitCellCoord(0, -3, 1, -2),
      xtal::UnitCellCoord(0, -3, 1, 4),
      xtal::UnitCellCoord(0, -3, 4, -2),
      xtal::UnitCellCoord(0, -3, 4, 1),
      xtal::UnitCellCoord(0, -2, -3, 1),
      xtal::UnitCellCoord(0, -2, -3, 4),
      xtal::UnitCellCoord(0, -2, 1, -3),
      xtal::UnitCellCoord(0, -2, 1, 4),
      xtal::UnitCellCoord(0, -2, 4, -3),
      xtal::UnitCellCoord(0, -2, 4, 1),
      xtal::UnitCellCoord(0, -1, -4, 2),
      xtal::UnitCellCoord(0, -1, -4, 3),
      xtal::UnitCellCoord(0, -1, 2, -4),
      xtal::UnitCellCoord(0, -1, 2, 3),
      xtal::UnitCellCoord(0, -1, 3, -4),
      xtal::UnitCellCoord(0, -1, 3, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -3, -2),
      xtal::UnitCellCoord(0, 1, -3, 4),
      xtal::UnitCellCoord(0, 1, -2, -3),
      xtal::UnitCellCoord(0, 1, -2, 4),
      xtal::UnitCellCoord(0, 1, 4, -3),
      xtal::UnitCellCoord(0, 1, 4, -2),
      xtal::UnitCellCoord(0, 2, -4, -1),
      xtal::UnitCellCoord(0, 2, -4, 3),
      xtal::UnitCellCoord(0, 2, -1, -4),
      xtal::UnitCellCoord(0, 2, -1, 3),
      xtal::UnitCellCoord(0, 2, 3, -4),
      xtal::UnitCellCoord(0, 2, 3, -1),
      xtal::UnitCellCoord(0, 3, -4, -1),
      xtal::UnitCellCoord(0, 3, -4, 2),
      xtal::UnitCellCoord(0, 3, -1, -4),
      xtal::UnitCellCoord(0, 3, -1, 2),
      xtal::UnitCellCoord(0, 3, 2, -4),
      xtal::UnitCellCoord(0, 3, 2, -1),
      xtal::UnitCellCoord(0, 4, -3, -2),
      xtal::UnitCellCoord(0, 4, -3, 1),
      xtal::UnitCellCoord(0, 4, -2, -3),
      xtal::UnitCellCoord(0, 4, -2, 1),
      xtal::UnitCellCoord(0, 4, 1, -3),
      xtal::UnitCellCoord(0, 4, 1, -2)
    };

    m_orbit_neighborhood[18] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-4, 0, 0),
      xtal::UnitCell(-4, 0, 4),
      xtal::UnitCell(-4, 4, 0),
      xtal::UnitCell(0, -4, 0),
      xtal::UnitCell(0, -4, 4),
      xtal::UnitCell(0, 0, -4),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 4),
      xtal::UnitCell(0, 4, -4),
      xtal::UnitCell(0, 4, 0),
      xtal::UnitCell(4, -4, 0),
      xtal::UnitCell(4, 0, -4),
      xtal::UnitCell(4, 0, 0)
    };

    m_orbit_site_neighborhood[18] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -4, 0, 0),
      xtal::UnitCellCoord(0, -4, 0, 4),
      xtal::UnitCellCoord(0, -4, 4, 0),
      xtal::UnitCellCoord(0, 0, -4, 0),
      xtal::UnitCellCoord(0, 0, -4, 4),
      xtal::UnitCellCoord(0, 0, 0, -4),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 4),
      xtal::UnitCellCoord(0, 0, 4, -4),
      xtal::UnitCellCoord(0, 0, 4, 0),
      xtal::UnitCellCoord(0, 4, -4, 0),
      xtal::UnitCellCoord(0, 4, 0, -4),
      xtal::UnitCellCoord(0, 4, 0, 0)
    };

    m_orbit_neighborhood[19] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-4, -1, 1),
      xtal::UnitCell(-4, -1, 4),
      xtal::UnitCell(-4, 1, -1),
      xtal::UnitCell(-4, 1, 4),
      xtal::UnitCell(-4, 4, -1),
      xtal::UnitCell(-4, 4, 1),
      xtal::UnitCell(-1, -4, 1),
      xtal::UnitCell(-1, -4, 4),
      xtal::UnitCell(-1, 1, -4),
      xtal::UnitCell(-1, 1, 4),
      xtal::UnitCell(-1, 4, -4),
      xtal::UnitCell(-1, 4, 1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -4, -1),
      xtal::UnitCell(1, -4, 4),
      xtal::UnitCell(1, -1, -4),
      xtal::UnitCell(1, -1, 4),
      xtal::UnitCell(1, 4, -4),
      xtal::UnitCell(1, 4, -1),
      xtal::UnitCell(4, -4, -1),
      xtal::UnitCell(4, -4, 1),
      xtal::UnitCell(4, -1, -4),
      xtal::UnitCell(4, -1, 1),
      xtal::UnitCell(4, 1, -4),
      xtal::UnitCell(4, 1, -1)
    };

    m_orbit_site_neighborhood[19] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -4, -1, 1),
      xtal::UnitCellCoord(0, -4, -1, 4),
      xtal::UnitCellCoord(0, -4, 1, -1),
      xtal::UnitCellCoord(0, -4, 1, 4),
      xtal::UnitCellCoord(0, -4, 4, -1),
      xtal::UnitCellCoord(0, -4, 4, 1),
      xtal::UnitCellCoord(0, -1, -4, 1),
      xtal::UnitCellCoord(0, -1, -4, 4),
      xtal::UnitCellCoord(0, -1, 1, -4),
      xtal::UnitCellCoord(0, -1, 1, 4),
      xtal::UnitCellCoord(0, -1, 4, -4),
      xtal::UnitCellCoord(0, -1, 4, 1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -4, -1),
      xtal::UnitCellCoord(0, 1, -4, 4),
      xtal::UnitCellCoord(0, 1, -1, -4),
      xtal::UnitCellCoord(0, 1, -1, 4),
      xtal::UnitCellCoord(0, 1, 4, -4),
      xtal::UnitCellCoord(0, 1, 4, -1),
      xtal::UnitCellCoord(0, 4, -4, -1),
      xtal::UnitCellCoord(0, 4, -4, 1),
      xtal::UnitCellCoord(0, 4, -1, -4),
      xtal::UnitCellCoord(0, 4, -1, 1),
      xtal::UnitCellCoord(0, 4, 1, -4),
      xtal::UnitCellCoord(0, 4, 1, -1)
    };

    m_orbit_neighborhood[20] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-5, 1, 2),
      xtal::UnitCell(-5, 2, 1),
      xtal::UnitCell(-5, 2, 2),
      xtal::UnitCell(-2, -2, -1),
      xtal::UnitCell(-2, -2, 5),
      xtal::UnitCell(-2, -1, -2),
      xtal::UnitCell(-2, -1, 5),
      xtal::UnitCell(-2, 5, -2),
      xtal::UnitCell(-2, 5, -1),
      xtal::UnitCell(-1, -2, -2),
      xtal::UnitCell(-1, -2, 5),
      xtal::UnitCell(-1, 5, -2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -5, 2),
      xtal::UnitCell(1, 2, -5),
      xtal::UnitCell(1, 2, 2),
      xtal::UnitCell(2, -5, 1),
      xtal::UnitCell(2, -5, 2),
      xtal::UnitCell(2, 1, -5),
      xtal::UnitCell(2, 1, 2),
      xtal::UnitCell(2, 2, -5),
      xtal::UnitCell(2, 2, 1),
      xtal::UnitCell(5, -2, -2),
      xtal::UnitCell(5, -2, -1),
      xtal::UnitCell(5, -1, -2)
    };

    m_orbit_site_neighborhood[20] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -5, 1, 2),
      xtal::UnitCellCoord(0, -5, 2, 1),
      xtal::UnitCellCoord(0, -5, 2, 2),
      xtal::UnitCellCoord(0, -2, -2, -1),
      xtal::UnitCellCoord(0, -2, -2, 5),
      xtal::UnitCellCoord(0, -2, -1, -2),
      xtal::UnitCellCoord(0, -2, -1, 5),
      xtal::UnitCellCoord(0, -2, 5, -2),
      xtal::UnitCellCoord(0, -2, 5, -1),
      xtal::UnitCellCoord(0, -1, -2, -2),
      xtal::UnitCellCoord(0, -1, -2, 5),
      xtal::UnitCellCoord(0, -1, 5, -2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -5, 2),
      xtal::UnitCellCoord(0, 1, 2, -5),
      xtal::UnitCellCoord(0, 1, 2, 2),
      xtal::UnitCellCoord(0, 2, -5, 1),
      xtal::UnitCellCoord(0, 2, -5, 2),
      xtal::UnitCellCoord(0, 2, 1, -5),
      xtal::UnitCellCoord(0, 2, 1, 2),
      xtal::UnitCellCoord(0, 2, 2, -5),
      xtal::UnitCellCoord(0, 2, 2, 1),
      xtal::UnitCellCoord(0, 5, -2, -2),
      xtal::UnitCellCoord(0, 5, -2, -1),
      xtal::UnitCellCoord(0, 5, -1, -2)
    };

    m_orbit_neighborhood[21] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, -3, 3),
      xtal::UnitCell(-3, 3, -3),
      xtal::UnitCell(-3, 3, 3),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(3, -3, -3),
      xtal::UnitCell(3, -3, 3),
      xtal::UnitCell(3, 3, -3)
    };

    m_orbit_site_neighborhood[21] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, -3, 3),
      xtal::UnitCellCoord(0, -3, 3, -3),
      xtal::UnitCellCoord(0, -3, 3, 3),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 3, -3, -3),
      xtal::UnitCellCoord(0, 3, -3, 3),
      xtal::UnitCellCoord(0, 3, 3, -3)
    };

    m_orbit_neighborhood[22] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-5, 1, 1),
      xtal::UnitCell(-5, 1, 3),
      xtal::UnitCell(-5, 3, 1),
      xtal::UnitCell(-3, -1, -1),
      xtal::UnitCell(-3, -1, 5),
      xtal::UnitCell(-3, 5, -1),
      xtal::UnitCell(-1, -3, -1),
      xtal::UnitCell(-1, -3, 5),
      xtal::UnitCell(-1, -1, -3),
      xtal::UnitCell(-1, -1, 5),
      xtal::UnitCell(-1, 5, -3),
      xtal::UnitCell(-1, 5, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -5, 1),
      xtal::UnitCell(1, -5, 3),
      xtal::UnitCell(1, 1, -5),
      xtal::UnitCell(1, 1, 3),
      xtal::UnitCell(1, 3, -5),
      xtal::UnitCell(1, 3, 1),
      xtal::UnitCell(3, -5, 1),
      xtal::UnitCell(3, 1, -5),
      xtal::UnitCell(3, 1, 1),
      xtal::UnitCell(5, -3, -1),
      xtal::UnitCell(5, -1, -3),
      xtal::UnitCell(5, -1, -1)
    };

    m_orbit_site_neighborhood[22] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -5, 1, 1),
      xtal::UnitCellCoord(0, -5, 1, 3),
      xtal::UnitCellCoord(0, -5, 3, 1),
      xtal::UnitCellCoord(0, -3, -1, -1),
      xtal::UnitCellCoord(0, -3, -1, 5),
      xtal::UnitCellCoord(0, -3, 5, -1),
      xtal::UnitCellCoord(0, -1, -3, -1),
      xtal::UnitCellCoord(0, -1, -3, 5),
      xtal::UnitCellCoord(0, -1, -1, -3),
      xtal::UnitCellCoord(0, -1, -1, 5),
      xtal::UnitCellCoord(0, -1, 5, -3),
      xtal::UnitCellCoord(0, -1, 5, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -5, 1),
      xtal::UnitCellCoord(0, 1, -5, 3),
      xtal::UnitCellCoord(0, 1, 1, -5),
      xtal::UnitCellCoord(0, 1, 1, 3),
      xtal::UnitCellCoord(0, 1, 3, -5),
      xtal::UnitCellCoord(0, 1, 3, 1),
      xtal::UnitCellCoord(0, 3, -5, 1),
      xtal::UnitCellCoord(0, 3, 1, -5),
      xtal::UnitCellCoord(0, 3, 1, 1),
      xtal::UnitCellCoord(0, 5, -3, -1),
      xtal::UnitCellCoord(0, 5, -1, -3),
      xtal::UnitCellCoord(0, 5, -1, -1)
    };

    m_orbit_neighborhood[23] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0)
    };

    m_orbit_site_neighborhood[23] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0)
    };

    m_orbit_neighborhood[24] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, -1)
    };

    m_orbit_site_neighborhood[24] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, -1)
    };

    m_orbit_neighborhood[25] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -1)
    };

    m_orbit_site_neighborhood[25] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -1)
    };

    m_orbit_neighborhood[26] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, -1),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -1)
    };

    m_orbit_site_neighborhood[26] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, -1),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -1)
    };

    m_orbit_neighborhood[27] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -1)
    };

    m_orbit_site_neighborhood[27] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -1)
    };

    m_orbit_neighborhood[28] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, -1),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -1)
    };

    m_orbit_site_neighborhood[28] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, -1),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -1)
    };

    m_orbit_neighborhood[29] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -1)
    };

    m_orbit_site_neighborhood[29] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -1)
    };

    m_orbit_neighborhood[30] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[30] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 2),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -2, 2),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -2),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 0, 2),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, -2),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 0, -2),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[31] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 1, -1),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[31] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 2),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -2, 2),
      xtal::UnitCellCoord(0, 0, 0, -2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 2),
      xtal::UnitCellCoord(0, 0, 2, -2),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 1, -1),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 0, -2),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[32] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, -1),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[32] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 0, 2),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -2, 2),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, -2),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 0, 2),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -2),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -2),
      xtal::UnitCellCoord(0, 2, 0, -1),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[33] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 1),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 1, 1),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, -1, 2),
      xtal::UnitCell(-1, 0, -1),
      xtal::UnitCell(-1, 0, 2),
      xtal::UnitCell(-1, 2, -1),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 1),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, -1, -1),
      xtal::UnitCell(0, -1, 2),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 1, -2),
      xtal::UnitCell(0, 1, 1),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, -1),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -2, 1),
      xtal::UnitCell(1, 0, -2),
      xtal::UnitCell(1, 0, 1),
      xtal::UnitCell(1, 1, -2),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, -1, -1),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, -1),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[33] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 1),
      xtal::UnitCellCoord(0, -2, 0, 2),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 1, 1),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, -1, 2),
      xtal::UnitCellCoord(0, -1, 0, -1),
      xtal::UnitCellCoord(0, -1, 0, 2),
      xtal::UnitCellCoord(0, -1, 2, -1),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -2, 1),
      xtal::UnitCellCoord(0, 0, -2, 2),
      xtal::UnitCellCoord(0, 0, -1, -1),
      xtal::UnitCellCoord(0, 0, -1, 2),
      xtal::UnitCellCoord(0, 0, 0, -2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 2),
      xtal::UnitCellCoord(0, 0, 1, -2),
      xtal::UnitCellCoord(0, 0, 1, 1),
      xtal::UnitCellCoord(0, 0, 2, -2),
      xtal::UnitCellCoord(0, 0, 2, -1),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -2, 1),
      xtal::UnitCellCoord(0, 1, 0, -2),
      xtal::UnitCellCoord(0, 1, 0, 1),
      xtal::UnitCellCoord(0, 1, 1, -2),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, -1, -1),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, -2),
      xtal::UnitCellCoord(0, 2, 0, -1),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[34] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 0, 2),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -2, 2),
      xtal::UnitCell(0, 0, -2),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 2),
      xtal::UnitCell(0, 2, -2),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 0, -2),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[34] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 2),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -2, 2),
      xtal::UnitCellCoord(0, 0, 0, -2),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 2),
      xtal::UnitCellCoord(0, 0, 2, -2),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 0, -2),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[35] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0)
    };

    m_orbit_site_neighborhood[35] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0)
    };

    m_orbit_neighborhood[36] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, -1)
    };

    m_orbit_site_neighborhood[36] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, -1)
    };

    m_orbit_neighborhood[37] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 1),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 0, 1),
      xtal::UnitCell(-1, 1, -1),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 1, 1),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, -1, 1),
      xtal::UnitCell(0, 0, -1),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 0, 1),
      xtal::UnitCell(0, 1, -1),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, -1),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, -1, 1),
      xtal::UnitCell(1, 0, -1),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, -1)
    };

    m_orbit_site_neighborhood[37] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 1),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 1),
      xtal::UnitCellCoord(0, -1, 1, -1),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 1, 1),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, -1, 1),
      xtal::UnitCellCoord(0, 0, 0, -1),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 0, 1),
      xtal::UnitCellCoord(0, 0, 1, -1),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, -1),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, -1, 1),
      xtal::UnitCellCoord(0, 1, 0, -1),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, -1)
    };

  }


  fcc_binary_Clexulator_default::~fcc_binary_Clexulator_default() {
    //nothing here for now
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void fcc_binary_Clexulator_default::_calc_global_corr_contribution(double *corr_begin) const {
    _calc_global_corr_contribution();
    for(size_type i = 0; i < corr_size(); i++) {
      *(corr_begin + i) = ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
    }
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void fcc_binary_Clexulator_default::_calc_global_corr_contribution() const {
    m_params.pre_eval();
    {
      _global_prepare<double>();
      for(size_type i = 0; i < corr_size(); i++) {
        ParamPack::Val<double>::set(m_params, m_corr_param_key, i, (this->*m_orbit_func_table_0[i])());
      }
    }
    m_params.post_eval();
  }

  /// \brief Calculate contribution to select global correlations from one unit cell
  void fcc_binary_Clexulator_default::_calc_restricted_global_corr_contribution(double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    _calc_restricted_global_corr_contribution(ind_list_begin, ind_list_end);
    for(; ind_list_begin < ind_list_end; ind_list_begin++) {
      *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(m_params, m_corr_param_key, *ind_list_begin);
    }
  }

  /// \brief Calculate contribution to select global correlations from one unit cell
  void fcc_binary_Clexulator_default::_calc_restricted_global_corr_contribution(size_type const *ind_list_begin, size_type const *ind_list_end) const {
    m_params.pre_eval();
    {
      _global_prepare<double>();
      for(; ind_list_begin < ind_list_end; ind_list_begin++) {
        ParamPack::Val<double>::set(m_params, m_corr_param_key, *ind_list_begin, (this->*m_orbit_func_table_0[*ind_list_begin])());
      }
    }
    m_params.post_eval();
  }

  /// \brief Calculate point correlations about basis site 'nlist_ind'
  void fcc_binary_Clexulator_default::_calc_point_corr(int nlist_ind, double *corr_begin) const {
    _calc_point_corr(nlist_ind);
    for(size_type i = 0; i < corr_size(); i++) {
      *(corr_begin + i) = ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
    }
  }

  /// \brief Calculate point correlations about basis site 'nlist_ind'
  void fcc_binary_Clexulator_default::_calc_point_corr(int nlist_ind) const {
    m_params.pre_eval();
    {
      _point_prepare<double>(nlist_ind);
      for(size_type i = 0; i < corr_size(); i++) {
        ParamPack::Val<double>::set(m_params, m_corr_param_key, i, (this->*m_flower_func_table_0[nlist_ind][i])());
      }
    }
    m_params.post_eval();
  }

  /// \brief Calculate select point correlations about basis site 'nlist_ind'
  void fcc_binary_Clexulator_default::_calc_restricted_point_corr(int nlist_ind, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    _calc_restricted_point_corr(nlist_ind, ind_list_begin, ind_list_end);
    for(; ind_list_begin < ind_list_end; ind_list_begin++) {
      *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(m_params, m_corr_param_key, *ind_list_begin);
    }
  }

  /// \brief Calculate select point correlations about basis site 'nlist_ind'
  void fcc_binary_Clexulator_default::_calc_restricted_point_corr(int nlist_ind, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    m_params.pre_eval();
    {
      _point_prepare<double>(nlist_ind);
      for(; ind_list_begin < ind_list_end; ind_list_begin++) {
        ParamPack::Val<double>::set(m_params, m_corr_param_key, *ind_list_begin, (this->*m_flower_func_table_0[nlist_ind][*ind_list_begin])());
      }
    }
    m_params.post_eval();
  }

  /// \brief Calculate the change in point correlations due to changing an occupant
  void fcc_binary_Clexulator_default::_calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double *corr_begin) const {
    _calc_delta_point_corr(nlist_ind, occ_i, occ_f);
    for(size_type i = 0; i < corr_size(); i++) {
      *(corr_begin + i) = ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
    }
  }

  /// \brief Calculate the change in point correlations due to changing an occupant
  void fcc_binary_Clexulator_default::_calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f) const {
    m_params.pre_eval();
    {
      _point_prepare<double>(nlist_ind);
     for(size_type i = 0; i < corr_size(); i++) {
        ParamPack::Val<double>::set(m_params, m_corr_param_key, i, (this->*m_delta_func_table_0[nlist_ind][i])(occ_i, occ_f));
      }
    }
    m_params.post_eval();
  }

  /// \brief Calculate the change in select point correlations due to changing an occupant
  void fcc_binary_Clexulator_default::_calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    _calc_restricted_delta_point_corr(nlist_ind, occ_i, occ_f, ind_list_begin, ind_list_end);
    for(; ind_list_begin < ind_list_end; ind_list_begin++) {
      *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(m_params, m_corr_param_key, *ind_list_begin);
    }
  }

  /// \brief Calculate the change in select point correlations due to changing an occupant
  void fcc_binary_Clexulator_default::_calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    m_params.pre_eval();
    {
      _point_prepare<double>(nlist_ind);
      for(; ind_list_begin < ind_list_end; ind_list_begin++) {
        ParamPack::Val<double>::set(m_params, m_corr_param_key, *ind_list_begin, (this->*m_delta_func_table_0[nlist_ind][*ind_list_begin])(occ_i, occ_f));
      }
    }
    m_params.post_eval();
  }


  template<typename Scalar>
  void fcc_binary_Clexulator_default::_point_prepare(int nlist_ind) const {
    switch(nlist_ind) {
    case 0:
      if(m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 0, eval_occ_func_0_0(0));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 1, eval_occ_func_0_0(1));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 2, eval_occ_func_0_0(2));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 3, eval_occ_func_0_0(3));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 4, eval_occ_func_0_0(4));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 5, eval_occ_func_0_0(5));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 6, eval_occ_func_0_0(6));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 7, eval_occ_func_0_0(7));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 8, eval_occ_func_0_0(8));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 9, eval_occ_func_0_0(9));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 10, eval_occ_func_0_0(10));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 11, eval_occ_func_0_0(11));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 12, eval_occ_func_0_0(12));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 13, eval_occ_func_0_0(13));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 14, eval_occ_func_0_0(14));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 15, eval_occ_func_0_0(15));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 16, eval_occ_func_0_0(16));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 17, eval_occ_func_0_0(17));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 18, eval_occ_func_0_0(18));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 19, eval_occ_func_0_0(19));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 20, eval_occ_func_0_0(20));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 21, eval_occ_func_0_0(21));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 22, eval_occ_func_0_0(22));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 23, eval_occ_func_0_0(23));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 24, eval_occ_func_0_0(24));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 25, eval_occ_func_0_0(25));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 26, eval_occ_func_0_0(26));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 27, eval_occ_func_0_0(27));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 28, eval_occ_func_0_0(28));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 29, eval_occ_func_0_0(29));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 30, eval_occ_func_0_0(30));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 31, eval_occ_func_0_0(31));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 32, eval_occ_func_0_0(32));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 33, eval_occ_func_0_0(33));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 34, eval_occ_func_0_0(34));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 35, eval_occ_func_0_0(35));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 36, eval_occ_func_0_0(36));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 37, eval_occ_func_0_0(37));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 38, eval_occ_func_0_0(38));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 39, eval_occ_func_0_0(39));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 40, eval_occ_func_0_0(40));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 41, eval_occ_func_0_0(41));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 42, eval_occ_func_0_0(42));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 43, eval_occ_func_0_0(43));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 44, eval_occ_func_0_0(44));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 45, eval_occ_func_0_0(45));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 46, eval_occ_func_0_0(46));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 47, eval_occ_func_0_0(47));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 48, eval_occ_func_0_0(48));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 49, eval_occ_func_0_0(49));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 50, eval_occ_func_0_0(50));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 51, eval_occ_func_0_0(51));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 52, eval_occ_func_0_0(52));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 53, eval_occ_func_0_0(53));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 54, eval_occ_func_0_0(54));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 55, eval_occ_func_0_0(55));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 56, eval_occ_func_0_0(56));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 57, eval_occ_func_0_0(57));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 58, eval_occ_func_0_0(58));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 59, eval_occ_func_0_0(59));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 60, eval_occ_func_0_0(60));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 61, eval_occ_func_0_0(61));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 62, eval_occ_func_0_0(62));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 63, eval_occ_func_0_0(63));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 64, eval_occ_func_0_0(64));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 65, eval_occ_func_0_0(65));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 66, eval_occ_func_0_0(66));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 67, eval_occ_func_0_0(67));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 68, eval_occ_func_0_0(68));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 69, eval_occ_func_0_0(69));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 70, eval_occ_func_0_0(70));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 71, eval_occ_func_0_0(71));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 72, eval_occ_func_0_0(72));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 73, eval_occ_func_0_0(73));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 74, eval_occ_func_0_0(74));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 75, eval_occ_func_0_0(75));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 76, eval_occ_func_0_0(76));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 77, eval_occ_func_0_0(77));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 78, eval_occ_func_0_0(78));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 79, eval_occ_func_0_0(79));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 80, eval_occ_func_0_0(80));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 81, eval_occ_func_0_0(81));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 82, eval_occ_func_0_0(82));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 83, eval_occ_func_0_0(83));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 84, eval_occ_func_0_0(84));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 85, eval_occ_func_0_0(85));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 86, eval_occ_func_0_0(86));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 87, eval_occ_func_0_0(87));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 88, eval_occ_func_0_0(88));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 89, eval_occ_func_0_0(89));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 90, eval_occ_func_0_0(90));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 91, eval_occ_func_0_0(91));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 92, eval_occ_func_0_0(92));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 93, eval_occ_func_0_0(93));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 94, eval_occ_func_0_0(94));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 95, eval_occ_func_0_0(95));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 96, eval_occ_func_0_0(96));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 97, eval_occ_func_0_0(97));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 98, eval_occ_func_0_0(98));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 99, eval_occ_func_0_0(99));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 100, eval_occ_func_0_0(100));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 101, eval_occ_func_0_0(101));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 102, eval_occ_func_0_0(102));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 103, eval_occ_func_0_0(103));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 104, eval_occ_func_0_0(104));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 105, eval_occ_func_0_0(105));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 106, eval_occ_func_0_0(106));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 107, eval_occ_func_0_0(107));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 108, eval_occ_func_0_0(108));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 109, eval_occ_func_0_0(109));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 110, eval_occ_func_0_0(110));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 111, eval_occ_func_0_0(111));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 112, eval_occ_func_0_0(112));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 113, eval_occ_func_0_0(113));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 114, eval_occ_func_0_0(114));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 115, eval_occ_func_0_0(115));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 116, eval_occ_func_0_0(116));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 117, eval_occ_func_0_0(117));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 118, eval_occ_func_0_0(118));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 119, eval_occ_func_0_0(119));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 120, eval_occ_func_0_0(120));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 121, eval_occ_func_0_0(121));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 122, eval_occ_func_0_0(122));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 123, eval_occ_func_0_0(123));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 124, eval_occ_func_0_0(124));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 125, eval_occ_func_0_0(125));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 126, eval_occ_func_0_0(126));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 127, eval_occ_func_0_0(127));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 128, eval_occ_func_0_0(128));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 129, eval_occ_func_0_0(129));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 130, eval_occ_func_0_0(130));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 131, eval_occ_func_0_0(131));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 132, eval_occ_func_0_0(132));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 133, eval_occ_func_0_0(133));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 134, eval_occ_func_0_0(134));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 135, eval_occ_func_0_0(135));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 136, eval_occ_func_0_0(136));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 137, eval_occ_func_0_0(137));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 138, eval_occ_func_0_0(138));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 139, eval_occ_func_0_0(139));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 140, eval_occ_func_0_0(140));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 141, eval_occ_func_0_0(141));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 142, eval_occ_func_0_0(142));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 143, eval_occ_func_0_0(143));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 144, eval_occ_func_0_0(144));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 145, eval_occ_func_0_0(145));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 146, eval_occ_func_0_0(146));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 147, eval_occ_func_0_0(147));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 148, eval_occ_func_0_0(148));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 149, eval_occ_func_0_0(149));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 150, eval_occ_func_0_0(150));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 151, eval_occ_func_0_0(151));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 152, eval_occ_func_0_0(152));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 153, eval_occ_func_0_0(153));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 154, eval_occ_func_0_0(154));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 155, eval_occ_func_0_0(155));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 156, eval_occ_func_0_0(156));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 157, eval_occ_func_0_0(157));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 158, eval_occ_func_0_0(158));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 159, eval_occ_func_0_0(159));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 160, eval_occ_func_0_0(160));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 161, eval_occ_func_0_0(161));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 162, eval_occ_func_0_0(162));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 163, eval_occ_func_0_0(163));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 164, eval_occ_func_0_0(164));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 165, eval_occ_func_0_0(165));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 166, eval_occ_func_0_0(166));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 167, eval_occ_func_0_0(167));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 168, eval_occ_func_0_0(168));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 169, eval_occ_func_0_0(169));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 170, eval_occ_func_0_0(170));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 171, eval_occ_func_0_0(171));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 172, eval_occ_func_0_0(172));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 173, eval_occ_func_0_0(173));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 174, eval_occ_func_0_0(174));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 175, eval_occ_func_0_0(175));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 176, eval_occ_func_0_0(176));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 177, eval_occ_func_0_0(177));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 178, eval_occ_func_0_0(178));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 179, eval_occ_func_0_0(179));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 180, eval_occ_func_0_0(180));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 181, eval_occ_func_0_0(181));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 182, eval_occ_func_0_0(182));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 183, eval_occ_func_0_0(183));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 184, eval_occ_func_0_0(184));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 185, eval_occ_func_0_0(185));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 186, eval_occ_func_0_0(186));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 187, eval_occ_func_0_0(187));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 188, eval_occ_func_0_0(188));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 189, eval_occ_func_0_0(189));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 190, eval_occ_func_0_0(190));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 191, eval_occ_func_0_0(191));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 192, eval_occ_func_0_0(192));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 193, eval_occ_func_0_0(193));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 194, eval_occ_func_0_0(194));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 195, eval_occ_func_0_0(195));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 196, eval_occ_func_0_0(196));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 197, eval_occ_func_0_0(197));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 198, eval_occ_func_0_0(198));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 199, eval_occ_func_0_0(199));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 200, eval_occ_func_0_0(200));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 201, eval_occ_func_0_0(201));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 202, eval_occ_func_0_0(202));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 203, eval_occ_func_0_0(203));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 204, eval_occ_func_0_0(204));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 205, eval_occ_func_0_0(205));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 206, eval_occ_func_0_0(206));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 207, eval_occ_func_0_0(207));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 208, eval_occ_func_0_0(208));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 209, eval_occ_func_0_0(209));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 210, eval_occ_func_0_0(210));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 211, eval_occ_func_0_0(211));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 212, eval_occ_func_0_0(212));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 213, eval_occ_func_0_0(213));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 214, eval_occ_func_0_0(214));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 215, eval_occ_func_0_0(215));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 216, eval_occ_func_0_0(216));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 217, eval_occ_func_0_0(217));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 218, eval_occ_func_0_0(218));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 219, eval_occ_func_0_0(219));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 220, eval_occ_func_0_0(220));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 221, eval_occ_func_0_0(221));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 222, eval_occ_func_0_0(222));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 223, eval_occ_func_0_0(223));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 224, eval_occ_func_0_0(224));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 225, eval_occ_func_0_0(225));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 226, eval_occ_func_0_0(226));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 227, eval_occ_func_0_0(227));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 228, eval_occ_func_0_0(228));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 229, eval_occ_func_0_0(229));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 230, eval_occ_func_0_0(230));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 231, eval_occ_func_0_0(231));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 232, eval_occ_func_0_0(232));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 233, eval_occ_func_0_0(233));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 234, eval_occ_func_0_0(234));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 235, eval_occ_func_0_0(235));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 236, eval_occ_func_0_0(236));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 237, eval_occ_func_0_0(237));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 238, eval_occ_func_0_0(238));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 239, eval_occ_func_0_0(239));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 240, eval_occ_func_0_0(240));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 241, eval_occ_func_0_0(241));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 242, eval_occ_func_0_0(242));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 243, eval_occ_func_0_0(243));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 244, eval_occ_func_0_0(244));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 245, eval_occ_func_0_0(245));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 246, eval_occ_func_0_0(246));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 247, eval_occ_func_0_0(247));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 248, eval_occ_func_0_0(248));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 249, eval_occ_func_0_0(249));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 250, eval_occ_func_0_0(250));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 251, eval_occ_func_0_0(251));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 252, eval_occ_func_0_0(252));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 253, eval_occ_func_0_0(253));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 254, eval_occ_func_0_0(254));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 255, eval_occ_func_0_0(255));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 256, eval_occ_func_0_0(256));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 257, eval_occ_func_0_0(257));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 258, eval_occ_func_0_0(258));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 259, eval_occ_func_0_0(259));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 260, eval_occ_func_0_0(260));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 261, eval_occ_func_0_0(261));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 262, eval_occ_func_0_0(262));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 263, eval_occ_func_0_0(263));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 264, eval_occ_func_0_0(264));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 265, eval_occ_func_0_0(265));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 266, eval_occ_func_0_0(266));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 267, eval_occ_func_0_0(267));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 268, eval_occ_func_0_0(268));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 269, eval_occ_func_0_0(269));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 270, eval_occ_func_0_0(270));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 271, eval_occ_func_0_0(271));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 272, eval_occ_func_0_0(272));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 273, eval_occ_func_0_0(273));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 274, eval_occ_func_0_0(274));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 275, eval_occ_func_0_0(275));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 276, eval_occ_func_0_0(276));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 277, eval_occ_func_0_0(277));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 278, eval_occ_func_0_0(278));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 279, eval_occ_func_0_0(279));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 280, eval_occ_func_0_0(280));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 281, eval_occ_func_0_0(281));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 282, eval_occ_func_0_0(282));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 283, eval_occ_func_0_0(283));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 284, eval_occ_func_0_0(284));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 285, eval_occ_func_0_0(285));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 286, eval_occ_func_0_0(286));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 287, eval_occ_func_0_0(287));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 288, eval_occ_func_0_0(288));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 289, eval_occ_func_0_0(289));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 290, eval_occ_func_0_0(290));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 291, eval_occ_func_0_0(291));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 292, eval_occ_func_0_0(292));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 293, eval_occ_func_0_0(293));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 294, eval_occ_func_0_0(294));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 295, eval_occ_func_0_0(295));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 296, eval_occ_func_0_0(296));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 297, eval_occ_func_0_0(297));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 298, eval_occ_func_0_0(298));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 299, eval_occ_func_0_0(299));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 300, eval_occ_func_0_0(300));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 301, eval_occ_func_0_0(301));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 302, eval_occ_func_0_0(302));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 303, eval_occ_func_0_0(303));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 304, eval_occ_func_0_0(304));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 305, eval_occ_func_0_0(305));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 306, eval_occ_func_0_0(306));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 307, eval_occ_func_0_0(307));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 308, eval_occ_func_0_0(308));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 309, eval_occ_func_0_0(309));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 310, eval_occ_func_0_0(310));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 311, eval_occ_func_0_0(311));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 312, eval_occ_func_0_0(312));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 313, eval_occ_func_0_0(313));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 314, eval_occ_func_0_0(314));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 315, eval_occ_func_0_0(315));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 316, eval_occ_func_0_0(316));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 317, eval_occ_func_0_0(317));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 318, eval_occ_func_0_0(318));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 319, eval_occ_func_0_0(319));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 320, eval_occ_func_0_0(320));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 321, eval_occ_func_0_0(321));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 322, eval_occ_func_0_0(322));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 323, eval_occ_func_0_0(323));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 324, eval_occ_func_0_0(324));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 325, eval_occ_func_0_0(325));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 326, eval_occ_func_0_0(326));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 327, eval_occ_func_0_0(327));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 328, eval_occ_func_0_0(328));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 329, eval_occ_func_0_0(329));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 330, eval_occ_func_0_0(330));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 331, eval_occ_func_0_0(331));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 332, eval_occ_func_0_0(332));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 333, eval_occ_func_0_0(333));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 334, eval_occ_func_0_0(334));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 335, eval_occ_func_0_0(335));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 336, eval_occ_func_0_0(336));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 337, eval_occ_func_0_0(337));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 338, eval_occ_func_0_0(338));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 339, eval_occ_func_0_0(339));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 340, eval_occ_func_0_0(340));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 341, eval_occ_func_0_0(341));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 342, eval_occ_func_0_0(342));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 343, eval_occ_func_0_0(343));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 344, eval_occ_func_0_0(344));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 345, eval_occ_func_0_0(345));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 346, eval_occ_func_0_0(346));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 347, eval_occ_func_0_0(347));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 348, eval_occ_func_0_0(348));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 349, eval_occ_func_0_0(349));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 350, eval_occ_func_0_0(350));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 351, eval_occ_func_0_0(351));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 352, eval_occ_func_0_0(352));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 353, eval_occ_func_0_0(353));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 354, eval_occ_func_0_0(354));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 355, eval_occ_func_0_0(355));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 356, eval_occ_func_0_0(356));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 357, eval_occ_func_0_0(357));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 358, eval_occ_func_0_0(358));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 359, eval_occ_func_0_0(359));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 360, eval_occ_func_0_0(360));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 361, eval_occ_func_0_0(361));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 362, eval_occ_func_0_0(362));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 363, eval_occ_func_0_0(363));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 364, eval_occ_func_0_0(364));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 365, eval_occ_func_0_0(365));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 366, eval_occ_func_0_0(366));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 367, eval_occ_func_0_0(367));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 368, eval_occ_func_0_0(368));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 369, eval_occ_func_0_0(369));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 370, eval_occ_func_0_0(370));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 371, eval_occ_func_0_0(371));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 372, eval_occ_func_0_0(372));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 373, eval_occ_func_0_0(373));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 374, eval_occ_func_0_0(374));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 375, eval_occ_func_0_0(375));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 376, eval_occ_func_0_0(376));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 377, eval_occ_func_0_0(377));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 378, eval_occ_func_0_0(378));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 379, eval_occ_func_0_0(379));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 380, eval_occ_func_0_0(380));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 381, eval_occ_func_0_0(381));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 382, eval_occ_func_0_0(382));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 383, eval_occ_func_0_0(383));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 384, eval_occ_func_0_0(384));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 385, eval_occ_func_0_0(385));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 386, eval_occ_func_0_0(386));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 387, eval_occ_func_0_0(387));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 388, eval_occ_func_0_0(388));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 389, eval_occ_func_0_0(389));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 390, eval_occ_func_0_0(390));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 391, eval_occ_func_0_0(391));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 392, eval_occ_func_0_0(392));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 393, eval_occ_func_0_0(393));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 394, eval_occ_func_0_0(394));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 395, eval_occ_func_0_0(395));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 396, eval_occ_func_0_0(396));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 397, eval_occ_func_0_0(397));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 398, eval_occ_func_0_0(398));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 399, eval_occ_func_0_0(399));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 400, eval_occ_func_0_0(400));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 401, eval_occ_func_0_0(401));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 402, eval_occ_func_0_0(402));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 403, eval_occ_func_0_0(403));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 404, eval_occ_func_0_0(404));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 405, eval_occ_func_0_0(405));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 406, eval_occ_func_0_0(406));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 407, eval_occ_func_0_0(407));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 408, eval_occ_func_0_0(408));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 409, eval_occ_func_0_0(409));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 410, eval_occ_func_0_0(410));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 411, eval_occ_func_0_0(411));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 412, eval_occ_func_0_0(412));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 413, eval_occ_func_0_0(413));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 414, eval_occ_func_0_0(414));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 415, eval_occ_func_0_0(415));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 416, eval_occ_func_0_0(416));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 417, eval_occ_func_0_0(417));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 418, eval_occ_func_0_0(418));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 419, eval_occ_func_0_0(419));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 420, eval_occ_func_0_0(420));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 421, eval_occ_func_0_0(421));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 422, eval_occ_func_0_0(422));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 423, eval_occ_func_0_0(423));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 424, eval_occ_func_0_0(424));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 425, eval_occ_func_0_0(425));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 426, eval_occ_func_0_0(426));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 427, eval_occ_func_0_0(427));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 428, eval_occ_func_0_0(428));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 429, eval_occ_func_0_0(429));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 430, eval_occ_func_0_0(430));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 431, eval_occ_func_0_0(431));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 432, eval_occ_func_0_0(432));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 433, eval_occ_func_0_0(433));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 434, eval_occ_func_0_0(434));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 435, eval_occ_func_0_0(435));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 436, eval_occ_func_0_0(436));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 437, eval_occ_func_0_0(437));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 438, eval_occ_func_0_0(438));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 439, eval_occ_func_0_0(439));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 440, eval_occ_func_0_0(440));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 441, eval_occ_func_0_0(441));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 442, eval_occ_func_0_0(442));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 443, eval_occ_func_0_0(443));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 444, eval_occ_func_0_0(444));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 445, eval_occ_func_0_0(445));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 446, eval_occ_func_0_0(446));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 447, eval_occ_func_0_0(447));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 448, eval_occ_func_0_0(448));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 449, eval_occ_func_0_0(449));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 450, eval_occ_func_0_0(450));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 451, eval_occ_func_0_0(451));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 452, eval_occ_func_0_0(452));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 453, eval_occ_func_0_0(453));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 454, eval_occ_func_0_0(454));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 455, eval_occ_func_0_0(455));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 456, eval_occ_func_0_0(456));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 457, eval_occ_func_0_0(457));
        ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 458, eval_occ_func_0_0(458));
      }
      break;
    }
  }
  template<typename Scalar>
  void fcc_binary_Clexulator_default::_global_prepare() const {
  if(m_params.eval_mode(m_occ_site_func_param_key) != ParamPack::READ) {
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 0, eval_occ_func_0_0(0));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 1, eval_occ_func_0_0(1));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 2, eval_occ_func_0_0(2));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 3, eval_occ_func_0_0(3));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 4, eval_occ_func_0_0(4));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 5, eval_occ_func_0_0(5));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 6, eval_occ_func_0_0(6));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 7, eval_occ_func_0_0(7));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 8, eval_occ_func_0_0(8));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 9, eval_occ_func_0_0(9));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 10, eval_occ_func_0_0(10));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 11, eval_occ_func_0_0(11));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 12, eval_occ_func_0_0(12));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 13, eval_occ_func_0_0(13));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 14, eval_occ_func_0_0(14));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 15, eval_occ_func_0_0(15));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 16, eval_occ_func_0_0(16));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 17, eval_occ_func_0_0(17));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 18, eval_occ_func_0_0(18));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 19, eval_occ_func_0_0(19));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 20, eval_occ_func_0_0(20));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 21, eval_occ_func_0_0(21));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 22, eval_occ_func_0_0(22));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 23, eval_occ_func_0_0(23));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 24, eval_occ_func_0_0(24));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 25, eval_occ_func_0_0(25));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 26, eval_occ_func_0_0(26));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 27, eval_occ_func_0_0(27));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 28, eval_occ_func_0_0(28));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 29, eval_occ_func_0_0(29));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 30, eval_occ_func_0_0(30));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 31, eval_occ_func_0_0(31));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 32, eval_occ_func_0_0(32));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 33, eval_occ_func_0_0(33));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 34, eval_occ_func_0_0(34));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 35, eval_occ_func_0_0(35));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 36, eval_occ_func_0_0(36));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 37, eval_occ_func_0_0(37));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 38, eval_occ_func_0_0(38));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 39, eval_occ_func_0_0(39));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 40, eval_occ_func_0_0(40));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 41, eval_occ_func_0_0(41));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 42, eval_occ_func_0_0(42));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 43, eval_occ_func_0_0(43));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 44, eval_occ_func_0_0(44));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 45, eval_occ_func_0_0(45));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 46, eval_occ_func_0_0(46));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 47, eval_occ_func_0_0(47));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 48, eval_occ_func_0_0(48));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 49, eval_occ_func_0_0(49));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 50, eval_occ_func_0_0(50));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 51, eval_occ_func_0_0(51));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 52, eval_occ_func_0_0(52));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 53, eval_occ_func_0_0(53));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 54, eval_occ_func_0_0(54));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 55, eval_occ_func_0_0(55));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 56, eval_occ_func_0_0(56));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 57, eval_occ_func_0_0(57));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 58, eval_occ_func_0_0(58));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 59, eval_occ_func_0_0(59));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 60, eval_occ_func_0_0(60));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 61, eval_occ_func_0_0(61));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 62, eval_occ_func_0_0(62));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 63, eval_occ_func_0_0(63));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 64, eval_occ_func_0_0(64));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 65, eval_occ_func_0_0(65));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 66, eval_occ_func_0_0(66));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 67, eval_occ_func_0_0(67));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 68, eval_occ_func_0_0(68));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 69, eval_occ_func_0_0(69));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 70, eval_occ_func_0_0(70));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 71, eval_occ_func_0_0(71));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 72, eval_occ_func_0_0(72));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 73, eval_occ_func_0_0(73));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 74, eval_occ_func_0_0(74));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 75, eval_occ_func_0_0(75));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 76, eval_occ_func_0_0(76));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 77, eval_occ_func_0_0(77));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 78, eval_occ_func_0_0(78));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 79, eval_occ_func_0_0(79));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 80, eval_occ_func_0_0(80));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 81, eval_occ_func_0_0(81));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 82, eval_occ_func_0_0(82));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 83, eval_occ_func_0_0(83));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 84, eval_occ_func_0_0(84));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 85, eval_occ_func_0_0(85));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 86, eval_occ_func_0_0(86));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 87, eval_occ_func_0_0(87));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 88, eval_occ_func_0_0(88));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 89, eval_occ_func_0_0(89));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 90, eval_occ_func_0_0(90));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 91, eval_occ_func_0_0(91));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 92, eval_occ_func_0_0(92));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 93, eval_occ_func_0_0(93));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 94, eval_occ_func_0_0(94));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 95, eval_occ_func_0_0(95));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 96, eval_occ_func_0_0(96));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 97, eval_occ_func_0_0(97));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 98, eval_occ_func_0_0(98));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 99, eval_occ_func_0_0(99));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 100, eval_occ_func_0_0(100));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 101, eval_occ_func_0_0(101));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 102, eval_occ_func_0_0(102));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 103, eval_occ_func_0_0(103));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 104, eval_occ_func_0_0(104));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 105, eval_occ_func_0_0(105));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 106, eval_occ_func_0_0(106));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 107, eval_occ_func_0_0(107));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 108, eval_occ_func_0_0(108));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 109, eval_occ_func_0_0(109));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 110, eval_occ_func_0_0(110));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 111, eval_occ_func_0_0(111));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 112, eval_occ_func_0_0(112));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 113, eval_occ_func_0_0(113));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 114, eval_occ_func_0_0(114));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 115, eval_occ_func_0_0(115));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 116, eval_occ_func_0_0(116));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 117, eval_occ_func_0_0(117));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 118, eval_occ_func_0_0(118));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 119, eval_occ_func_0_0(119));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 120, eval_occ_func_0_0(120));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 121, eval_occ_func_0_0(121));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 122, eval_occ_func_0_0(122));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 123, eval_occ_func_0_0(123));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 124, eval_occ_func_0_0(124));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 125, eval_occ_func_0_0(125));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 126, eval_occ_func_0_0(126));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 127, eval_occ_func_0_0(127));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 128, eval_occ_func_0_0(128));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 129, eval_occ_func_0_0(129));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 130, eval_occ_func_0_0(130));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 131, eval_occ_func_0_0(131));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 132, eval_occ_func_0_0(132));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 133, eval_occ_func_0_0(133));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 134, eval_occ_func_0_0(134));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 135, eval_occ_func_0_0(135));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 136, eval_occ_func_0_0(136));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 137, eval_occ_func_0_0(137));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 138, eval_occ_func_0_0(138));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 139, eval_occ_func_0_0(139));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 140, eval_occ_func_0_0(140));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 141, eval_occ_func_0_0(141));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 142, eval_occ_func_0_0(142));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 143, eval_occ_func_0_0(143));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 144, eval_occ_func_0_0(144));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 145, eval_occ_func_0_0(145));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 146, eval_occ_func_0_0(146));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 147, eval_occ_func_0_0(147));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 148, eval_occ_func_0_0(148));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 149, eval_occ_func_0_0(149));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 150, eval_occ_func_0_0(150));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 151, eval_occ_func_0_0(151));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 152, eval_occ_func_0_0(152));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 153, eval_occ_func_0_0(153));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 154, eval_occ_func_0_0(154));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 155, eval_occ_func_0_0(155));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 156, eval_occ_func_0_0(156));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 157, eval_occ_func_0_0(157));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 158, eval_occ_func_0_0(158));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 159, eval_occ_func_0_0(159));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 160, eval_occ_func_0_0(160));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 161, eval_occ_func_0_0(161));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 162, eval_occ_func_0_0(162));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 163, eval_occ_func_0_0(163));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 164, eval_occ_func_0_0(164));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 165, eval_occ_func_0_0(165));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 166, eval_occ_func_0_0(166));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 167, eval_occ_func_0_0(167));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 168, eval_occ_func_0_0(168));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 169, eval_occ_func_0_0(169));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 170, eval_occ_func_0_0(170));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 171, eval_occ_func_0_0(171));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 172, eval_occ_func_0_0(172));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 173, eval_occ_func_0_0(173));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 174, eval_occ_func_0_0(174));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 175, eval_occ_func_0_0(175));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 176, eval_occ_func_0_0(176));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 177, eval_occ_func_0_0(177));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 178, eval_occ_func_0_0(178));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 179, eval_occ_func_0_0(179));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 180, eval_occ_func_0_0(180));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 181, eval_occ_func_0_0(181));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 182, eval_occ_func_0_0(182));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 183, eval_occ_func_0_0(183));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 184, eval_occ_func_0_0(184));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 185, eval_occ_func_0_0(185));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 186, eval_occ_func_0_0(186));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 187, eval_occ_func_0_0(187));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 188, eval_occ_func_0_0(188));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 189, eval_occ_func_0_0(189));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 190, eval_occ_func_0_0(190));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 191, eval_occ_func_0_0(191));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 192, eval_occ_func_0_0(192));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 193, eval_occ_func_0_0(193));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 194, eval_occ_func_0_0(194));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 195, eval_occ_func_0_0(195));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 196, eval_occ_func_0_0(196));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 197, eval_occ_func_0_0(197));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 198, eval_occ_func_0_0(198));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 199, eval_occ_func_0_0(199));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 200, eval_occ_func_0_0(200));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 201, eval_occ_func_0_0(201));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 202, eval_occ_func_0_0(202));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 203, eval_occ_func_0_0(203));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 204, eval_occ_func_0_0(204));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 205, eval_occ_func_0_0(205));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 206, eval_occ_func_0_0(206));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 207, eval_occ_func_0_0(207));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 208, eval_occ_func_0_0(208));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 209, eval_occ_func_0_0(209));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 210, eval_occ_func_0_0(210));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 211, eval_occ_func_0_0(211));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 212, eval_occ_func_0_0(212));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 213, eval_occ_func_0_0(213));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 214, eval_occ_func_0_0(214));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 215, eval_occ_func_0_0(215));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 216, eval_occ_func_0_0(216));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 217, eval_occ_func_0_0(217));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 218, eval_occ_func_0_0(218));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 219, eval_occ_func_0_0(219));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 220, eval_occ_func_0_0(220));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 221, eval_occ_func_0_0(221));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 222, eval_occ_func_0_0(222));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 223, eval_occ_func_0_0(223));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 224, eval_occ_func_0_0(224));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 225, eval_occ_func_0_0(225));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 226, eval_occ_func_0_0(226));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 227, eval_occ_func_0_0(227));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 228, eval_occ_func_0_0(228));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 229, eval_occ_func_0_0(229));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 230, eval_occ_func_0_0(230));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 231, eval_occ_func_0_0(231));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 232, eval_occ_func_0_0(232));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 233, eval_occ_func_0_0(233));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 234, eval_occ_func_0_0(234));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 235, eval_occ_func_0_0(235));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 236, eval_occ_func_0_0(236));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 237, eval_occ_func_0_0(237));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 238, eval_occ_func_0_0(238));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 239, eval_occ_func_0_0(239));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 240, eval_occ_func_0_0(240));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 241, eval_occ_func_0_0(241));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 242, eval_occ_func_0_0(242));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 243, eval_occ_func_0_0(243));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 244, eval_occ_func_0_0(244));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 245, eval_occ_func_0_0(245));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 246, eval_occ_func_0_0(246));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 247, eval_occ_func_0_0(247));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 248, eval_occ_func_0_0(248));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 249, eval_occ_func_0_0(249));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 250, eval_occ_func_0_0(250));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 251, eval_occ_func_0_0(251));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 252, eval_occ_func_0_0(252));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 253, eval_occ_func_0_0(253));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 254, eval_occ_func_0_0(254));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 255, eval_occ_func_0_0(255));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 256, eval_occ_func_0_0(256));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 257, eval_occ_func_0_0(257));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 258, eval_occ_func_0_0(258));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 259, eval_occ_func_0_0(259));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 260, eval_occ_func_0_0(260));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 261, eval_occ_func_0_0(261));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 262, eval_occ_func_0_0(262));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 263, eval_occ_func_0_0(263));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 264, eval_occ_func_0_0(264));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 265, eval_occ_func_0_0(265));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 266, eval_occ_func_0_0(266));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 267, eval_occ_func_0_0(267));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 268, eval_occ_func_0_0(268));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 269, eval_occ_func_0_0(269));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 270, eval_occ_func_0_0(270));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 271, eval_occ_func_0_0(271));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 272, eval_occ_func_0_0(272));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 273, eval_occ_func_0_0(273));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 274, eval_occ_func_0_0(274));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 275, eval_occ_func_0_0(275));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 276, eval_occ_func_0_0(276));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 277, eval_occ_func_0_0(277));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 278, eval_occ_func_0_0(278));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 279, eval_occ_func_0_0(279));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 280, eval_occ_func_0_0(280));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 281, eval_occ_func_0_0(281));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 282, eval_occ_func_0_0(282));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 283, eval_occ_func_0_0(283));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 284, eval_occ_func_0_0(284));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 285, eval_occ_func_0_0(285));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 286, eval_occ_func_0_0(286));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 287, eval_occ_func_0_0(287));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 288, eval_occ_func_0_0(288));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 289, eval_occ_func_0_0(289));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 290, eval_occ_func_0_0(290));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 291, eval_occ_func_0_0(291));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 292, eval_occ_func_0_0(292));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 293, eval_occ_func_0_0(293));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 294, eval_occ_func_0_0(294));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 295, eval_occ_func_0_0(295));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 296, eval_occ_func_0_0(296));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 297, eval_occ_func_0_0(297));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 298, eval_occ_func_0_0(298));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 299, eval_occ_func_0_0(299));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 300, eval_occ_func_0_0(300));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 301, eval_occ_func_0_0(301));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 302, eval_occ_func_0_0(302));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 303, eval_occ_func_0_0(303));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 304, eval_occ_func_0_0(304));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 305, eval_occ_func_0_0(305));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 306, eval_occ_func_0_0(306));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 307, eval_occ_func_0_0(307));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 308, eval_occ_func_0_0(308));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 309, eval_occ_func_0_0(309));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 310, eval_occ_func_0_0(310));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 311, eval_occ_func_0_0(311));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 312, eval_occ_func_0_0(312));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 313, eval_occ_func_0_0(313));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 314, eval_occ_func_0_0(314));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 315, eval_occ_func_0_0(315));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 316, eval_occ_func_0_0(316));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 317, eval_occ_func_0_0(317));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 318, eval_occ_func_0_0(318));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 319, eval_occ_func_0_0(319));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 320, eval_occ_func_0_0(320));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 321, eval_occ_func_0_0(321));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 322, eval_occ_func_0_0(322));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 323, eval_occ_func_0_0(323));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 324, eval_occ_func_0_0(324));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 325, eval_occ_func_0_0(325));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 326, eval_occ_func_0_0(326));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 327, eval_occ_func_0_0(327));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 328, eval_occ_func_0_0(328));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 329, eval_occ_func_0_0(329));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 330, eval_occ_func_0_0(330));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 331, eval_occ_func_0_0(331));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 332, eval_occ_func_0_0(332));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 333, eval_occ_func_0_0(333));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 334, eval_occ_func_0_0(334));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 335, eval_occ_func_0_0(335));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 336, eval_occ_func_0_0(336));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 337, eval_occ_func_0_0(337));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 338, eval_occ_func_0_0(338));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 339, eval_occ_func_0_0(339));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 340, eval_occ_func_0_0(340));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 341, eval_occ_func_0_0(341));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 342, eval_occ_func_0_0(342));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 343, eval_occ_func_0_0(343));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 344, eval_occ_func_0_0(344));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 345, eval_occ_func_0_0(345));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 346, eval_occ_func_0_0(346));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 347, eval_occ_func_0_0(347));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 348, eval_occ_func_0_0(348));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 349, eval_occ_func_0_0(349));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 350, eval_occ_func_0_0(350));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 351, eval_occ_func_0_0(351));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 352, eval_occ_func_0_0(352));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 353, eval_occ_func_0_0(353));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 354, eval_occ_func_0_0(354));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 355, eval_occ_func_0_0(355));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 356, eval_occ_func_0_0(356));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 357, eval_occ_func_0_0(357));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 358, eval_occ_func_0_0(358));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 359, eval_occ_func_0_0(359));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 360, eval_occ_func_0_0(360));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 361, eval_occ_func_0_0(361));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 362, eval_occ_func_0_0(362));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 363, eval_occ_func_0_0(363));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 364, eval_occ_func_0_0(364));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 365, eval_occ_func_0_0(365));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 366, eval_occ_func_0_0(366));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 367, eval_occ_func_0_0(367));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 368, eval_occ_func_0_0(368));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 369, eval_occ_func_0_0(369));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 370, eval_occ_func_0_0(370));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 371, eval_occ_func_0_0(371));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 372, eval_occ_func_0_0(372));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 373, eval_occ_func_0_0(373));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 374, eval_occ_func_0_0(374));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 375, eval_occ_func_0_0(375));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 376, eval_occ_func_0_0(376));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 377, eval_occ_func_0_0(377));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 378, eval_occ_func_0_0(378));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 379, eval_occ_func_0_0(379));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 380, eval_occ_func_0_0(380));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 381, eval_occ_func_0_0(381));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 382, eval_occ_func_0_0(382));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 383, eval_occ_func_0_0(383));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 384, eval_occ_func_0_0(384));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 385, eval_occ_func_0_0(385));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 386, eval_occ_func_0_0(386));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 387, eval_occ_func_0_0(387));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 388, eval_occ_func_0_0(388));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 389, eval_occ_func_0_0(389));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 390, eval_occ_func_0_0(390));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 391, eval_occ_func_0_0(391));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 392, eval_occ_func_0_0(392));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 393, eval_occ_func_0_0(393));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 394, eval_occ_func_0_0(394));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 395, eval_occ_func_0_0(395));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 396, eval_occ_func_0_0(396));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 397, eval_occ_func_0_0(397));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 398, eval_occ_func_0_0(398));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 399, eval_occ_func_0_0(399));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 400, eval_occ_func_0_0(400));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 401, eval_occ_func_0_0(401));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 402, eval_occ_func_0_0(402));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 403, eval_occ_func_0_0(403));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 404, eval_occ_func_0_0(404));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 405, eval_occ_func_0_0(405));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 406, eval_occ_func_0_0(406));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 407, eval_occ_func_0_0(407));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 408, eval_occ_func_0_0(408));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 409, eval_occ_func_0_0(409));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 410, eval_occ_func_0_0(410));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 411, eval_occ_func_0_0(411));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 412, eval_occ_func_0_0(412));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 413, eval_occ_func_0_0(413));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 414, eval_occ_func_0_0(414));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 415, eval_occ_func_0_0(415));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 416, eval_occ_func_0_0(416));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 417, eval_occ_func_0_0(417));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 418, eval_occ_func_0_0(418));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 419, eval_occ_func_0_0(419));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 420, eval_occ_func_0_0(420));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 421, eval_occ_func_0_0(421));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 422, eval_occ_func_0_0(422));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 423, eval_occ_func_0_0(423));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 424, eval_occ_func_0_0(424));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 425, eval_occ_func_0_0(425));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 426, eval_occ_func_0_0(426));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 427, eval_occ_func_0_0(427));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 428, eval_occ_func_0_0(428));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 429, eval_occ_func_0_0(429));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 430, eval_occ_func_0_0(430));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 431, eval_occ_func_0_0(431));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 432, eval_occ_func_0_0(432));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 433, eval_occ_func_0_0(433));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 434, eval_occ_func_0_0(434));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 435, eval_occ_func_0_0(435));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 436, eval_occ_func_0_0(436));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 437, eval_occ_func_0_0(437));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 438, eval_occ_func_0_0(438));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 439, eval_occ_func_0_0(439));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 440, eval_occ_func_0_0(440));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 441, eval_occ_func_0_0(441));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 442, eval_occ_func_0_0(442));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 443, eval_occ_func_0_0(443));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 444, eval_occ_func_0_0(444));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 445, eval_occ_func_0_0(445));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 446, eval_occ_func_0_0(446));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 447, eval_occ_func_0_0(447));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 448, eval_occ_func_0_0(448));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 449, eval_occ_func_0_0(449));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 450, eval_occ_func_0_0(450));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 451, eval_occ_func_0_0(451));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 452, eval_occ_func_0_0(452));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 453, eval_occ_func_0_0(453));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 454, eval_occ_func_0_0(454));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 455, eval_occ_func_0_0(455));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 456, eval_occ_func_0_0(456));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 457, eval_occ_func_0_0(457));
    ParamPack::Val<Scalar>::set(m_params, m_occ_site_func_param_key, 0, 458, eval_occ_func_0_0(458));
  }
  }

  // Basis functions for empty cluster:
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_0_0() const {
    return 1;
  }

  /**** Basis functions for orbit 1****
0.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_1_0() const {
    return occ_func_0_0(0);
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_1_0_at_0() const {
    return occ_func_0_0(0);
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_1_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]);
  }

  /**** Basis functions for orbit 2****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_2_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(0) * occ_func_0_0(8)) / 6.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_2_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(0)) / 6.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_2_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(9) + occ_func_0_0(4) + occ_func_0_0(10) + occ_func_0_0(3) + occ_func_0_0(11) + occ_func_0_0(2) + occ_func_0_0(12) + occ_func_0_0(1) + occ_func_0_0(7) + occ_func_0_0(6) + occ_func_0_0(8) + occ_func_0_0(5)) / 6.;
  }

  /**** Basis functions for orbit 3****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -1.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_3_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(17) * occ_func_0_0(0)) / 3.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_3_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(14)) / 3.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_3_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(16) + occ_func_0_0(15) + occ_func_0_0(18) + occ_func_0_0(13) + occ_func_0_0(17) + occ_func_0_0(14)) / 3.;
  }

  /**** Basis functions for orbit 4****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_4_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(35) + occ_func_0_0(0) * occ_func_0_0(42) + occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(32) * occ_func_0_0(0) + occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(31) * occ_func_0_0(0) + occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(34) + occ_func_0_0(0) * occ_func_0_0(38) + occ_func_0_0(33) * occ_func_0_0(0) + occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(37) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_4_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(42) + occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(32) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(29) + occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(22) + occ_func_0_0(31) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(30) + occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25) + occ_func_0_0(0) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(38) + occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(33) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(28) + occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(37) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(24)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_4_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(35) + occ_func_0_0(26) + occ_func_0_0(42) + occ_func_0_0(19) + occ_func_0_0(41) + occ_func_0_0(20) + occ_func_0_0(32) + occ_func_0_0(29) + occ_func_0_0(39) + occ_func_0_0(22) + occ_func_0_0(31) + occ_func_0_0(30) + occ_func_0_0(36) + occ_func_0_0(25) + occ_func_0_0(34) + occ_func_0_0(27) + occ_func_0_0(38) + occ_func_0_0(23) + occ_func_0_0(33) + occ_func_0_0(28) + occ_func_0_0(40) + occ_func_0_0(21) + occ_func_0_0(37) + occ_func_0_0(24)) / 12.;
  }

  /**** Basis functions for orbit 5****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_5_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(0) * occ_func_0_0(50)) / 6.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_5_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(46) * occ_func_0_0(0) + occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(44) + occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(43) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(48) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(47) * occ_func_0_0(0)) / 6.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_5_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(51) + occ_func_0_0(46) + occ_func_0_0(52) + occ_func_0_0(45) + occ_func_0_0(53) + occ_func_0_0(44) + occ_func_0_0(54) + occ_func_0_0(43) + occ_func_0_0(49) + occ_func_0_0(48) + occ_func_0_0(50) + occ_func_0_0(47)) / 6.;
  }

  /**** Basis functions for orbit 6****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 2.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_6_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(68) + occ_func_0_0(0) * occ_func_0_0(78) + occ_func_0_0(0) * occ_func_0_0(74) + occ_func_0_0(0) * occ_func_0_0(67) + occ_func_0_0(0) * occ_func_0_0(70) + occ_func_0_0(75) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(72) + occ_func_0_0(69) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(73) + occ_func_0_0(0) * occ_func_0_0(77) + occ_func_0_0(71) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(76)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_6_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(68) + occ_func_0_0(65) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(78) + occ_func_0_0(55) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(74) + occ_func_0_0(59) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(67) + occ_func_0_0(66) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(70) + occ_func_0_0(63) * occ_func_0_0(0) + occ_func_0_0(75) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(58) + occ_func_0_0(0) * occ_func_0_0(72) + occ_func_0_0(61) * occ_func_0_0(0) + occ_func_0_0(69) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(64) + occ_func_0_0(0) * occ_func_0_0(73) + occ_func_0_0(60) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(77) + occ_func_0_0(56) * occ_func_0_0(0) + occ_func_0_0(71) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(62) + occ_func_0_0(0) * occ_func_0_0(76) + occ_func_0_0(57) * occ_func_0_0(0)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_6_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(68) + occ_func_0_0(65) + occ_func_0_0(78) + occ_func_0_0(55) + occ_func_0_0(74) + occ_func_0_0(59) + occ_func_0_0(67) + occ_func_0_0(66) + occ_func_0_0(70) + occ_func_0_0(63) + occ_func_0_0(75) + occ_func_0_0(58) + occ_func_0_0(72) + occ_func_0_0(61) + occ_func_0_0(69) + occ_func_0_0(64) + occ_func_0_0(73) + occ_func_0_0(60) + occ_func_0_0(77) + occ_func_0_0(56) + occ_func_0_0(71) + occ_func_0_0(62) + occ_func_0_0(76) + occ_func_0_0(57)) / 12.;
  }

  /**** Basis functions for orbit 7****
0.0000000 0.0000000 0.0000000 A  B  

3.0000000 -1.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_7_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(86) + occ_func_0_0(0) * occ_func_0_0(84) + occ_func_0_0(0) * occ_func_0_0(85) + occ_func_0_0(0) * occ_func_0_0(83)) / 4.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_7_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(86) + occ_func_0_0(79) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(84) + occ_func_0_0(81) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(85) + occ_func_0_0(80) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(83) + occ_func_0_0(82) * occ_func_0_0(0)) / 4.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_7_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(86) + occ_func_0_0(79) + occ_func_0_0(84) + occ_func_0_0(81) + occ_func_0_0(85) + occ_func_0_0(80) + occ_func_0_0(83) + occ_func_0_0(82)) / 4.;
  }

  /**** Basis functions for orbit 8****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 -3.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_8_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(111) + occ_func_0_0(122) * occ_func_0_0(0) + occ_func_0_0(118) * occ_func_0_0(0) + occ_func_0_0(123) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(133) + occ_func_0_0(126) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(129) + occ_func_0_0(115) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(125) + occ_func_0_0(131) * occ_func_0_0(0) + occ_func_0_0(134) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(114) + occ_func_0_0(0) * occ_func_0_0(121) + occ_func_0_0(0) * occ_func_0_0(128) + occ_func_0_0(0) * occ_func_0_0(117) + occ_func_0_0(0) * occ_func_0_0(132) + occ_func_0_0(119) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(116) + occ_func_0_0(112) * occ_func_0_0(0) + occ_func_0_0(113) * occ_func_0_0(0) + occ_func_0_0(127) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(124) + occ_func_0_0(130) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(120)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_8_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(111) + occ_func_0_0(110) * occ_func_0_0(0) + occ_func_0_0(122) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(99) + occ_func_0_0(118) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(103) + occ_func_0_0(123) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(98) + occ_func_0_0(0) * occ_func_0_0(133) + occ_func_0_0(88) * occ_func_0_0(0) + occ_func_0_0(126) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(95) + occ_func_0_0(0) * occ_func_0_0(129) + occ_func_0_0(92) * occ_func_0_0(0) + occ_func_0_0(115) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(106) + occ_func_0_0(0) * occ_func_0_0(125) + occ_func_0_0(96) * occ_func_0_0(0) + occ_func_0_0(131) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(90) + occ_func_0_0(134) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(87) + occ_func_0_0(0) * occ_func_0_0(114) + occ_func_0_0(107) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(121) + occ_func_0_0(100) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(128) + occ_func_0_0(93) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(117) + occ_func_0_0(104) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(132) + occ_func_0_0(89) * occ_func_0_0(0) + occ_func_0_0(119) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(102) + occ_func_0_0(0) * occ_func_0_0(116) + occ_func_0_0(105) * occ_func_0_0(0) + occ_func_0_0(112) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(109) + occ_func_0_0(113) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(108) + occ_func_0_0(127) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(94) + occ_func_0_0(0) * occ_func_0_0(124) + occ_func_0_0(97) * occ_func_0_0(0) + occ_func_0_0(130) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(91) + occ_func_0_0(0) * occ_func_0_0(120) + occ_func_0_0(101) * occ_func_0_0(0)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_8_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(111) + occ_func_0_0(110) + occ_func_0_0(122) + occ_func_0_0(99) + occ_func_0_0(118) + occ_func_0_0(103) + occ_func_0_0(123) + occ_func_0_0(98) + occ_func_0_0(133) + occ_func_0_0(88) + occ_func_0_0(126) + occ_func_0_0(95) + occ_func_0_0(129) + occ_func_0_0(92) + occ_func_0_0(115) + occ_func_0_0(106) + occ_func_0_0(125) + occ_func_0_0(96) + occ_func_0_0(131) + occ_func_0_0(90) + occ_func_0_0(134) + occ_func_0_0(87) + occ_func_0_0(114) + occ_func_0_0(107) + occ_func_0_0(121) + occ_func_0_0(100) + occ_func_0_0(128) + occ_func_0_0(93) + occ_func_0_0(117) + occ_func_0_0(104) + occ_func_0_0(132) + occ_func_0_0(89) + occ_func_0_0(119) + occ_func_0_0(102) + occ_func_0_0(116) + occ_func_0_0(105) + occ_func_0_0(112) + occ_func_0_0(109) + occ_func_0_0(113) + occ_func_0_0(108) + occ_func_0_0(127) + occ_func_0_0(94) + occ_func_0_0(124) + occ_func_0_0(97) + occ_func_0_0(130) + occ_func_0_0(91) + occ_func_0_0(120) + occ_func_0_0(101)) / 24.;
  }

  /**** Basis functions for orbit 9****
0.0000000 0.0000000 0.0000000 A  B  

2.0000000 -2.0000000 -2.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_9_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(138) + occ_func_0_0(0) * occ_func_0_0(140) + occ_func_0_0(139) * occ_func_0_0(0)) / 3.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_9_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(138) + occ_func_0_0(137) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(140) + occ_func_0_0(135) * occ_func_0_0(0) + occ_func_0_0(139) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(136)) / 3.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_9_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(138) + occ_func_0_0(137) + occ_func_0_0(140) + occ_func_0_0(135) + occ_func_0_0(139) + occ_func_0_0(136)) / 3.;
  }

  /**** Basis functions for orbit 10****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 3.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_10_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(161) + occ_func_0_0(171) * occ_func_0_0(0) + occ_func_0_0(174) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(175) + occ_func_0_0(0) * occ_func_0_0(159) + occ_func_0_0(0) * occ_func_0_0(160)) / 6.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_10_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(161) + occ_func_0_0(156) * occ_func_0_0(0) + occ_func_0_0(171) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(146) + occ_func_0_0(174) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(143) + occ_func_0_0(0) * occ_func_0_0(175) + occ_func_0_0(142) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(159) + occ_func_0_0(158) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(160) + occ_func_0_0(157) * occ_func_0_0(0)) / 6.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_10_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(161) + occ_func_0_0(156) + occ_func_0_0(171) + occ_func_0_0(146) + occ_func_0_0(174) + occ_func_0_0(143) + occ_func_0_0(175) + occ_func_0_0(142) + occ_func_0_0(159) + occ_func_0_0(158) + occ_func_0_0(160) + occ_func_0_0(157)) / 6.;
  }

  /**** Basis functions for orbit 11****
0.0000000 0.0000000 0.0000000 A  B  

2.0000000 -3.0000000 2.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_11_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(166) + occ_func_0_0(0) * occ_func_0_0(176) + occ_func_0_0(0) * occ_func_0_0(173) + occ_func_0_0(0) * occ_func_0_0(162) + occ_func_0_0(170) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(163) + occ_func_0_0(167) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(165) + occ_func_0_0(0) * occ_func_0_0(169) + occ_func_0_0(164) * occ_func_0_0(0) + occ_func_0_0(172) * occ_func_0_0(0) + occ_func_0_0(168) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_11_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(166) + occ_func_0_0(151) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(176) + occ_func_0_0(141) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(173) + occ_func_0_0(144) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(162) + occ_func_0_0(155) * occ_func_0_0(0) + occ_func_0_0(170) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(147) + occ_func_0_0(0) * occ_func_0_0(163) + occ_func_0_0(154) * occ_func_0_0(0) + occ_func_0_0(167) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(150) + occ_func_0_0(0) * occ_func_0_0(165) + occ_func_0_0(152) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(169) + occ_func_0_0(148) * occ_func_0_0(0) + occ_func_0_0(164) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(153) + occ_func_0_0(172) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(145) + occ_func_0_0(168) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(149)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_11_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(166) + occ_func_0_0(151) + occ_func_0_0(176) + occ_func_0_0(141) + occ_func_0_0(173) + occ_func_0_0(144) + occ_func_0_0(162) + occ_func_0_0(155) + occ_func_0_0(170) + occ_func_0_0(147) + occ_func_0_0(163) + occ_func_0_0(154) + occ_func_0_0(167) + occ_func_0_0(150) + occ_func_0_0(165) + occ_func_0_0(152) + occ_func_0_0(169) + occ_func_0_0(148) + occ_func_0_0(164) + occ_func_0_0(153) + occ_func_0_0(172) + occ_func_0_0(145) + occ_func_0_0(168) + occ_func_0_0(149)) / 12.;
  }

  /**** Basis functions for orbit 12****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -3.0000000 3.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_12_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(190) + occ_func_0_0(0) * occ_func_0_0(200) + occ_func_0_0(0) * occ_func_0_0(196) + occ_func_0_0(0) * occ_func_0_0(189) + occ_func_0_0(0) * occ_func_0_0(192) + occ_func_0_0(197) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(194) + occ_func_0_0(191) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(195) + occ_func_0_0(0) * occ_func_0_0(199) + occ_func_0_0(193) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(198)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_12_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(190) + occ_func_0_0(187) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(200) + occ_func_0_0(177) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(196) + occ_func_0_0(181) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(189) + occ_func_0_0(188) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(192) + occ_func_0_0(185) * occ_func_0_0(0) + occ_func_0_0(197) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(180) + occ_func_0_0(0) * occ_func_0_0(194) + occ_func_0_0(183) * occ_func_0_0(0) + occ_func_0_0(191) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(186) + occ_func_0_0(0) * occ_func_0_0(195) + occ_func_0_0(182) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(199) + occ_func_0_0(178) * occ_func_0_0(0) + occ_func_0_0(193) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(184) + occ_func_0_0(0) * occ_func_0_0(198) + occ_func_0_0(179) * occ_func_0_0(0)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_12_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(190) + occ_func_0_0(187) + occ_func_0_0(200) + occ_func_0_0(177) + occ_func_0_0(196) + occ_func_0_0(181) + occ_func_0_0(189) + occ_func_0_0(188) + occ_func_0_0(192) + occ_func_0_0(185) + occ_func_0_0(197) + occ_func_0_0(180) + occ_func_0_0(194) + occ_func_0_0(183) + occ_func_0_0(191) + occ_func_0_0(186) + occ_func_0_0(195) + occ_func_0_0(182) + occ_func_0_0(199) + occ_func_0_0(178) + occ_func_0_0(193) + occ_func_0_0(184) + occ_func_0_0(198) + occ_func_0_0(179)) / 12.;
  }

  /**** Basis functions for orbit 13****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -4.0000000 1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_13_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(213) + occ_func_0_0(0) * occ_func_0_0(222) + occ_func_0_0(0) * occ_func_0_0(223) + occ_func_0_0(221) * occ_func_0_0(0) + occ_func_0_0(216) * occ_func_0_0(0) + occ_func_0_0(220) * occ_func_0_0(0) + occ_func_0_0(217) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(214) + occ_func_0_0(0) * occ_func_0_0(215) + occ_func_0_0(0) * occ_func_0_0(219) + occ_func_0_0(224) * occ_func_0_0(0) + occ_func_0_0(218) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_13_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(213) + occ_func_0_0(212) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(222) + occ_func_0_0(203) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(223) + occ_func_0_0(202) * occ_func_0_0(0) + occ_func_0_0(221) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(204) + occ_func_0_0(216) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(209) + occ_func_0_0(220) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(205) + occ_func_0_0(217) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(208) + occ_func_0_0(0) * occ_func_0_0(214) + occ_func_0_0(211) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(215) + occ_func_0_0(210) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(219) + occ_func_0_0(206) * occ_func_0_0(0) + occ_func_0_0(224) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(201) + occ_func_0_0(218) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(207)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_13_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(213) + occ_func_0_0(212) + occ_func_0_0(222) + occ_func_0_0(203) + occ_func_0_0(223) + occ_func_0_0(202) + occ_func_0_0(221) + occ_func_0_0(204) + occ_func_0_0(216) + occ_func_0_0(209) + occ_func_0_0(220) + occ_func_0_0(205) + occ_func_0_0(217) + occ_func_0_0(208) + occ_func_0_0(214) + occ_func_0_0(211) + occ_func_0_0(215) + occ_func_0_0(210) + occ_func_0_0(219) + occ_func_0_0(206) + occ_func_0_0(224) + occ_func_0_0(201) + occ_func_0_0(218) + occ_func_0_0(207)) / 12.;
  }

  /**** Basis functions for orbit 14****
0.0000000 0.0000000 0.0000000 A  B  

2.0000000 -4.0000000 2.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_14_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(241) + occ_func_0_0(0) * occ_func_0_0(248) + occ_func_0_0(0) * occ_func_0_0(247) + occ_func_0_0(238) * occ_func_0_0(0) + occ_func_0_0(245) * occ_func_0_0(0) + occ_func_0_0(237) * occ_func_0_0(0) + occ_func_0_0(242) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(240) + occ_func_0_0(0) * occ_func_0_0(244) + occ_func_0_0(239) * occ_func_0_0(0) + occ_func_0_0(246) * occ_func_0_0(0) + occ_func_0_0(243) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_14_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(241) + occ_func_0_0(232) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(248) + occ_func_0_0(225) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(247) + occ_func_0_0(226) * occ_func_0_0(0) + occ_func_0_0(238) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(235) + occ_func_0_0(245) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(228) + occ_func_0_0(237) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(236) + occ_func_0_0(242) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(231) + occ_func_0_0(0) * occ_func_0_0(240) + occ_func_0_0(233) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(244) + occ_func_0_0(229) * occ_func_0_0(0) + occ_func_0_0(239) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(234) + occ_func_0_0(246) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(227) + occ_func_0_0(243) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(230)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_14_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(241) + occ_func_0_0(232) + occ_func_0_0(248) + occ_func_0_0(225) + occ_func_0_0(247) + occ_func_0_0(226) + occ_func_0_0(238) + occ_func_0_0(235) + occ_func_0_0(245) + occ_func_0_0(228) + occ_func_0_0(237) + occ_func_0_0(236) + occ_func_0_0(242) + occ_func_0_0(231) + occ_func_0_0(240) + occ_func_0_0(233) + occ_func_0_0(244) + occ_func_0_0(229) + occ_func_0_0(239) + occ_func_0_0(234) + occ_func_0_0(246) + occ_func_0_0(227) + occ_func_0_0(243) + occ_func_0_0(230)) / 12.;
  }

  /**** Basis functions for orbit 15****
0.0000000 0.0000000 0.0000000 A  B  

2.0000000 -3.0000000 3.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_15_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(298) + occ_func_0_0(0) * occ_func_0_0(314) + occ_func_0_0(0) * occ_func_0_0(306) + occ_func_0_0(0) * occ_func_0_0(297) + occ_func_0_0(0) * occ_func_0_0(300) + occ_func_0_0(307) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(302) + occ_func_0_0(299) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(305) + occ_func_0_0(0) * occ_func_0_0(313) + occ_func_0_0(301) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(308)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_15_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(298) + occ_func_0_0(271) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(314) + occ_func_0_0(255) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(306) + occ_func_0_0(263) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(297) + occ_func_0_0(272) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(300) + occ_func_0_0(269) * occ_func_0_0(0) + occ_func_0_0(307) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(262) + occ_func_0_0(0) * occ_func_0_0(302) + occ_func_0_0(267) * occ_func_0_0(0) + occ_func_0_0(299) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(270) + occ_func_0_0(0) * occ_func_0_0(305) + occ_func_0_0(264) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(313) + occ_func_0_0(256) * occ_func_0_0(0) + occ_func_0_0(301) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(268) + occ_func_0_0(0) * occ_func_0_0(308) + occ_func_0_0(261) * occ_func_0_0(0)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_15_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(298) + occ_func_0_0(271) + occ_func_0_0(314) + occ_func_0_0(255) + occ_func_0_0(306) + occ_func_0_0(263) + occ_func_0_0(297) + occ_func_0_0(272) + occ_func_0_0(300) + occ_func_0_0(269) + occ_func_0_0(307) + occ_func_0_0(262) + occ_func_0_0(302) + occ_func_0_0(267) + occ_func_0_0(299) + occ_func_0_0(270) + occ_func_0_0(305) + occ_func_0_0(264) + occ_func_0_0(313) + occ_func_0_0(256) + occ_func_0_0(301) + occ_func_0_0(268) + occ_func_0_0(308) + occ_func_0_0(261)) / 12.;
  }

  /**** Basis functions for orbit 16****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 -4.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_16_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(285) + occ_func_0_0(296) * occ_func_0_0(0) + occ_func_0_0(292) * occ_func_0_0(0) + occ_func_0_0(303) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(319) + occ_func_0_0(310) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(315) + occ_func_0_0(289) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(309) + occ_func_0_0(317) * occ_func_0_0(0) + occ_func_0_0(320) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(288) + occ_func_0_0(0) * occ_func_0_0(295) + occ_func_0_0(0) * occ_func_0_0(312) + occ_func_0_0(0) * occ_func_0_0(291) + occ_func_0_0(0) * occ_func_0_0(318) + occ_func_0_0(293) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(290) + occ_func_0_0(286) * occ_func_0_0(0) + occ_func_0_0(287) * occ_func_0_0(0) + occ_func_0_0(311) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(304) + occ_func_0_0(316) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(294)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_16_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(285) + occ_func_0_0(284) * occ_func_0_0(0) + occ_func_0_0(296) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(273) + occ_func_0_0(292) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(277) + occ_func_0_0(303) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(266) + occ_func_0_0(0) * occ_func_0_0(319) + occ_func_0_0(250) * occ_func_0_0(0) + occ_func_0_0(310) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(259) + occ_func_0_0(0) * occ_func_0_0(315) + occ_func_0_0(254) * occ_func_0_0(0) + occ_func_0_0(289) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(280) + occ_func_0_0(0) * occ_func_0_0(309) + occ_func_0_0(260) * occ_func_0_0(0) + occ_func_0_0(317) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(252) + occ_func_0_0(320) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(249) + occ_func_0_0(0) * occ_func_0_0(288) + occ_func_0_0(281) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(295) + occ_func_0_0(274) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(312) + occ_func_0_0(257) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(291) + occ_func_0_0(278) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(318) + occ_func_0_0(251) * occ_func_0_0(0) + occ_func_0_0(293) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(276) + occ_func_0_0(0) * occ_func_0_0(290) + occ_func_0_0(279) * occ_func_0_0(0) + occ_func_0_0(286) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(283) + occ_func_0_0(287) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(282) + occ_func_0_0(311) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(258) + occ_func_0_0(0) * occ_func_0_0(304) + occ_func_0_0(265) * occ_func_0_0(0) + occ_func_0_0(316) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(253) + occ_func_0_0(0) * occ_func_0_0(294) + occ_func_0_0(275) * occ_func_0_0(0)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_16_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(285) + occ_func_0_0(284) + occ_func_0_0(296) + occ_func_0_0(273) + occ_func_0_0(292) + occ_func_0_0(277) + occ_func_0_0(303) + occ_func_0_0(266) + occ_func_0_0(319) + occ_func_0_0(250) + occ_func_0_0(310) + occ_func_0_0(259) + occ_func_0_0(315) + occ_func_0_0(254) + occ_func_0_0(289) + occ_func_0_0(280) + occ_func_0_0(309) + occ_func_0_0(260) + occ_func_0_0(317) + occ_func_0_0(252) + occ_func_0_0(320) + occ_func_0_0(249) + occ_func_0_0(288) + occ_func_0_0(281) + occ_func_0_0(295) + occ_func_0_0(274) + occ_func_0_0(312) + occ_func_0_0(257) + occ_func_0_0(291) + occ_func_0_0(278) + occ_func_0_0(318) + occ_func_0_0(251) + occ_func_0_0(293) + occ_func_0_0(276) + occ_func_0_0(290) + occ_func_0_0(279) + occ_func_0_0(286) + occ_func_0_0(283) + occ_func_0_0(287) + occ_func_0_0(282) + occ_func_0_0(311) + occ_func_0_0(258) + occ_func_0_0(304) + occ_func_0_0(265) + occ_func_0_0(316) + occ_func_0_0(253) + occ_func_0_0(294) + occ_func_0_0(275)) / 24.;
  }

  /**** Basis functions for orbit 17****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -3.0000000 -2.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_17_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(345) + occ_func_0_0(0) * occ_func_0_0(357) + occ_func_0_0(0) * occ_func_0_0(361) + occ_func_0_0(366) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(353) + occ_func_0_0(367) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(352) + occ_func_0_0(0) * occ_func_0_0(348) + occ_func_0_0(0) * occ_func_0_0(368) + occ_func_0_0(355) * occ_func_0_0(0) + occ_func_0_0(354) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(349) + occ_func_0_0(358) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(364) + occ_func_0_0(362) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(356) + occ_func_0_0(0) * occ_func_0_0(360) + occ_func_0_0(347) * occ_func_0_0(0) + occ_func_0_0(346) * occ_func_0_0(0) + occ_func_0_0(350) * occ_func_0_0(0) + occ_func_0_0(363) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(365) + occ_func_0_0(351) * occ_func_0_0(0) + occ_func_0_0(359) * occ_func_0_0(0)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_17_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(345) + occ_func_0_0(344) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(357) + occ_func_0_0(332) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(361) + occ_func_0_0(328) * occ_func_0_0(0) + occ_func_0_0(366) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(323) + occ_func_0_0(0) * occ_func_0_0(353) + occ_func_0_0(336) * occ_func_0_0(0) + occ_func_0_0(367) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(322) + occ_func_0_0(0) * occ_func_0_0(352) + occ_func_0_0(337) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(348) + occ_func_0_0(341) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(368) + occ_func_0_0(321) * occ_func_0_0(0) + occ_func_0_0(355) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(334) + occ_func_0_0(354) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(335) + occ_func_0_0(0) * occ_func_0_0(349) + occ_func_0_0(340) * occ_func_0_0(0) + occ_func_0_0(358) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(331) + occ_func_0_0(0) * occ_func_0_0(364) + occ_func_0_0(325) * occ_func_0_0(0) + occ_func_0_0(362) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(327) + occ_func_0_0(0) * occ_func_0_0(356) + occ_func_0_0(333) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(360) + occ_func_0_0(329) * occ_func_0_0(0) + occ_func_0_0(347) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(342) + occ_func_0_0(346) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(343) + occ_func_0_0(350) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(339) + occ_func_0_0(363) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(326) + occ_func_0_0(0) * occ_func_0_0(365) + occ_func_0_0(324) * occ_func_0_0(0) + occ_func_0_0(351) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(338) + occ_func_0_0(359) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(330)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_17_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(345) + occ_func_0_0(344) + occ_func_0_0(357) + occ_func_0_0(332) + occ_func_0_0(361) + occ_func_0_0(328) + occ_func_0_0(366) + occ_func_0_0(323) + occ_func_0_0(353) + occ_func_0_0(336) + occ_func_0_0(367) + occ_func_0_0(322) + occ_func_0_0(352) + occ_func_0_0(337) + occ_func_0_0(348) + occ_func_0_0(341) + occ_func_0_0(368) + occ_func_0_0(321) + occ_func_0_0(355) + occ_func_0_0(334) + occ_func_0_0(354) + occ_func_0_0(335) + occ_func_0_0(349) + occ_func_0_0(340) + occ_func_0_0(358) + occ_func_0_0(331) + occ_func_0_0(364) + occ_func_0_0(325) + occ_func_0_0(362) + occ_func_0_0(327) + occ_func_0_0(356) + occ_func_0_0(333) + occ_func_0_0(360) + occ_func_0_0(329) + occ_func_0_0(347) + occ_func_0_0(342) + occ_func_0_0(346) + occ_func_0_0(343) + occ_func_0_0(350) + occ_func_0_0(339) + occ_func_0_0(363) + occ_func_0_0(326) + occ_func_0_0(365) + occ_func_0_0(324) + occ_func_0_0(351) + occ_func_0_0(338) + occ_func_0_0(359) + occ_func_0_0(330)) / 24.;
  }

  /**** Basis functions for orbit 18****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 4.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_18_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(377) + occ_func_0_0(378) * occ_func_0_0(0) + occ_func_0_0(379) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(380) + occ_func_0_0(0) * occ_func_0_0(375) + occ_func_0_0(0) * occ_func_0_0(376)) / 6.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_18_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(377) + occ_func_0_0(372) * occ_func_0_0(0) + occ_func_0_0(378) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(371) + occ_func_0_0(379) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(370) + occ_func_0_0(0) * occ_func_0_0(380) + occ_func_0_0(369) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(375) + occ_func_0_0(374) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(376) + occ_func_0_0(373) * occ_func_0_0(0)) / 6.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_18_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(377) + occ_func_0_0(372) + occ_func_0_0(378) + occ_func_0_0(371) + occ_func_0_0(379) + occ_func_0_0(370) + occ_func_0_0(380) + occ_func_0_0(369) + occ_func_0_0(375) + occ_func_0_0(374) + occ_func_0_0(376) + occ_func_0_0(373)) / 6.;
  }

  /**** Basis functions for orbit 19****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -4.0000000 4.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_19_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(407) + occ_func_0_0(0) * occ_func_0_0(425) + occ_func_0_0(0) * occ_func_0_0(421) + occ_func_0_0(0) * occ_func_0_0(406) + occ_func_0_0(0) * occ_func_0_0(409) + occ_func_0_0(422) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(413) + occ_func_0_0(408) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(420) + occ_func_0_0(0) * occ_func_0_0(424) + occ_func_0_0(412) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(423)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_19_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(407) + occ_func_0_0(402) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(425) + occ_func_0_0(384) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(421) + occ_func_0_0(388) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(406) + occ_func_0_0(403) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(409) + occ_func_0_0(400) * occ_func_0_0(0) + occ_func_0_0(422) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(387) + occ_func_0_0(0) * occ_func_0_0(413) + occ_func_0_0(396) * occ_func_0_0(0) + occ_func_0_0(408) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(401) + occ_func_0_0(0) * occ_func_0_0(420) + occ_func_0_0(389) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(424) + occ_func_0_0(385) * occ_func_0_0(0) + occ_func_0_0(412) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(397) + occ_func_0_0(0) * occ_func_0_0(423) + occ_func_0_0(386) * occ_func_0_0(0)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_19_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(407) + occ_func_0_0(402) + occ_func_0_0(425) + occ_func_0_0(384) + occ_func_0_0(421) + occ_func_0_0(388) + occ_func_0_0(406) + occ_func_0_0(403) + occ_func_0_0(409) + occ_func_0_0(400) + occ_func_0_0(422) + occ_func_0_0(387) + occ_func_0_0(413) + occ_func_0_0(396) + occ_func_0_0(408) + occ_func_0_0(401) + occ_func_0_0(420) + occ_func_0_0(389) + occ_func_0_0(424) + occ_func_0_0(385) + occ_func_0_0(412) + occ_func_0_0(397) + occ_func_0_0(423) + occ_func_0_0(386)) / 12.;
  }

  /**** Basis functions for orbit 20****
0.0000000 0.0000000 0.0000000 A  B  

2.0000000 -5.0000000 2.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_20_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(415) + occ_func_0_0(0) * occ_func_0_0(428) + occ_func_0_0(0) * occ_func_0_0(427) + occ_func_0_0(411) * occ_func_0_0(0) + occ_func_0_0(419) * occ_func_0_0(0) + occ_func_0_0(410) * occ_func_0_0(0) + occ_func_0_0(416) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(414) + occ_func_0_0(0) * occ_func_0_0(418) + occ_func_0_0(0) * occ_func_0_0(405) + occ_func_0_0(426) * occ_func_0_0(0) + occ_func_0_0(417) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_20_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(415) + occ_func_0_0(394) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(428) + occ_func_0_0(381) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(427) + occ_func_0_0(382) * occ_func_0_0(0) + occ_func_0_0(411) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(398) + occ_func_0_0(419) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(390) + occ_func_0_0(410) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(399) + occ_func_0_0(416) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(393) + occ_func_0_0(0) * occ_func_0_0(414) + occ_func_0_0(395) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(418) + occ_func_0_0(391) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(405) + occ_func_0_0(404) * occ_func_0_0(0) + occ_func_0_0(426) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(383) + occ_func_0_0(417) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(392)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_20_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(415) + occ_func_0_0(394) + occ_func_0_0(428) + occ_func_0_0(381) + occ_func_0_0(427) + occ_func_0_0(382) + occ_func_0_0(411) + occ_func_0_0(398) + occ_func_0_0(419) + occ_func_0_0(390) + occ_func_0_0(410) + occ_func_0_0(399) + occ_func_0_0(416) + occ_func_0_0(393) + occ_func_0_0(414) + occ_func_0_0(395) + occ_func_0_0(418) + occ_func_0_0(391) + occ_func_0_0(405) + occ_func_0_0(404) + occ_func_0_0(426) + occ_func_0_0(383) + occ_func_0_0(417) + occ_func_0_0(392)) / 12.;
  }

  /**** Basis functions for orbit 21****
0.0000000 0.0000000 0.0000000 A  B  

3.0000000 -3.0000000 -3.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_21_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(451) + occ_func_0_0(0) * occ_func_0_0(455) + occ_func_0_0(452) * occ_func_0_0(0)) / 3.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_21_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(451) + occ_func_0_0(436) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(455) + occ_func_0_0(432) * occ_func_0_0(0) + occ_func_0_0(452) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(435)) / 3.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_21_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(451) + occ_func_0_0(436) + occ_func_0_0(455) + occ_func_0_0(432) + occ_func_0_0(452) + occ_func_0_0(435)) / 3.;
  }

  /**** Basis functions for orbit 22****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -5.0000000 1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_22_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(444) + occ_func_0_0(0) * occ_func_0_0(456) + occ_func_0_0(0) * occ_func_0_0(457) + occ_func_0_0(454) * occ_func_0_0(0) + occ_func_0_0(447) * occ_func_0_0(0) + occ_func_0_0(453) * occ_func_0_0(0) + occ_func_0_0(448) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(445) + occ_func_0_0(0) * occ_func_0_0(446) + occ_func_0_0(0) * occ_func_0_0(450) + occ_func_0_0(458) * occ_func_0_0(0) + occ_func_0_0(449) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_22_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(444) + occ_func_0_0(443) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(456) + occ_func_0_0(431) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(457) + occ_func_0_0(430) * occ_func_0_0(0) + occ_func_0_0(454) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(433) + occ_func_0_0(447) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(440) + occ_func_0_0(453) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(434) + occ_func_0_0(448) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(439) + occ_func_0_0(0) * occ_func_0_0(445) + occ_func_0_0(442) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(446) + occ_func_0_0(441) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(450) + occ_func_0_0(437) * occ_func_0_0(0) + occ_func_0_0(458) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(429) + occ_func_0_0(449) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(438)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_22_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(444) + occ_func_0_0(443) + occ_func_0_0(456) + occ_func_0_0(431) + occ_func_0_0(457) + occ_func_0_0(430) + occ_func_0_0(454) + occ_func_0_0(433) + occ_func_0_0(447) + occ_func_0_0(440) + occ_func_0_0(453) + occ_func_0_0(434) + occ_func_0_0(448) + occ_func_0_0(439) + occ_func_0_0(445) + occ_func_0_0(442) + occ_func_0_0(446) + occ_func_0_0(441) + occ_func_0_0(450) + occ_func_0_0(437) + occ_func_0_0(458) + occ_func_0_0(429) + occ_func_0_0(449) + occ_func_0_0(438)) / 12.;
  }

  /**** Basis functions for orbit 23****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 0.0000000 1.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_23_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(9) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(12) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(8) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(8)) / 8.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_23_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(9) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(3) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(11) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(2) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(10) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(3) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(7) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6)) / 8.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_23_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(9) + occ_func_0_0(6) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(5) + occ_func_0_0(10) * occ_func_0_0(12) + occ_func_0_0(9) * occ_func_0_0(3) + occ_func_0_0(4) * occ_func_0_0(1) + occ_func_0_0(11) * occ_func_0_0(10) + occ_func_0_0(8) * occ_func_0_0(3) + occ_func_0_0(5) * occ_func_0_0(2) + occ_func_0_0(11) * occ_func_0_0(12) + occ_func_0_0(2) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(1) + occ_func_0_0(12) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(6) + occ_func_0_0(2) * occ_func_0_0(1) + occ_func_0_0(11) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(2) * occ_func_0_0(3) + occ_func_0_0(12) * occ_func_0_0(9) + occ_func_0_0(10) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(3) + occ_func_0_0(9) * occ_func_0_0(8) + occ_func_0_0(7) * occ_func_0_0(5) + occ_func_0_0(4) * occ_func_0_0(6)) / 8.;
  }

  /**** Basis functions for orbit 24****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 0.0000000 1.0000000 A  B  

1.0000000 -1.0000000 1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_24_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(18) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(18) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(17)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_24_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(14) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(18) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(13) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(13) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(16) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(13) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(13) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(17) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(14) * occ_func_0_0(6) * occ_func_0_0(0)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_24_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(6) * occ_func_0_0(10) + occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(9) * occ_func_0_0(18) + occ_func_0_0(4) * occ_func_0_0(11) + occ_func_0_0(13) * occ_func_0_0(2) + occ_func_0_0(8) * occ_func_0_0(12) + occ_func_0_0(5) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(1) + occ_func_0_0(10) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(7) + occ_func_0_0(3) * occ_func_0_0(6) + occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(5) + occ_func_0_0(10) * occ_func_0_0(7) + occ_func_0_0(16) * occ_func_0_0(6) + occ_func_0_0(3) * occ_func_0_0(15) + occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(5) * occ_func_0_0(12) + occ_func_0_0(13) * occ_func_0_0(1) + occ_func_0_0(11) * occ_func_0_0(9) + occ_func_0_0(16) * occ_func_0_0(4) + occ_func_0_0(2) * occ_func_0_0(15) + occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(9) + occ_func_0_0(2) * occ_func_0_0(4) + occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(8) + occ_func_0_0(13) * occ_func_0_0(5) + occ_func_0_0(18) * occ_func_0_0(11) + occ_func_0_0(9) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(13) + occ_func_0_0(10) * occ_func_0_0(17) + occ_func_0_0(3) * occ_func_0_0(7) + occ_func_0_0(14) * occ_func_0_0(6)) / 12.;
  }

  /**** Basis functions for orbit 25****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 0.0000000 1.0000000 A  B  

1.0000000 0.0000000 1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_25_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(37) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(33) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(40) + occ_func_0_0(38) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(41) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(36) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(39) + occ_func_0_0(42) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(32) + occ_func_0_0(35) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(37) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(32) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(39) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(42) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(41) + occ_func_0_0(31) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(38) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(40) + occ_func_0_0(34) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(33)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_25_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(37) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(24) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(33) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(28) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(31) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(30) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(40) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(21) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(38) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(41) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(20) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(34) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(27) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(39) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(22) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(32) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(29) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(35) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(26) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(37) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(24) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(32) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(29) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(39) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(22) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(41) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(20) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(31) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(30) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(38) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(40) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(21) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(34) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(27) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(33) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(28) * occ_func_0_0(4) * occ_func_0_0(0)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_25_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(37) + occ_func_0_0(6) * occ_func_0_0(12) + occ_func_0_0(24) * occ_func_0_0(1) + occ_func_0_0(9) * occ_func_0_0(33) + occ_func_0_0(4) * occ_func_0_0(8) + occ_func_0_0(28) * occ_func_0_0(5) + occ_func_0_0(8) * occ_func_0_0(7) + occ_func_0_0(31) * occ_func_0_0(6) + occ_func_0_0(5) * occ_func_0_0(30) + occ_func_0_0(10) * occ_func_0_0(40) + occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(21) * occ_func_0_0(2) + occ_func_0_0(38) * occ_func_0_0(8) + occ_func_0_0(11) * occ_func_0_0(5) + occ_func_0_0(2) * occ_func_0_0(23) + occ_func_0_0(12) * occ_func_0_0(41) + occ_func_0_0(1) * occ_func_0_0(10) + occ_func_0_0(20) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(9) + occ_func_0_0(34) * occ_func_0_0(4) + occ_func_0_0(3) * occ_func_0_0(27) + occ_func_0_0(11) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(7) + occ_func_0_0(2) * occ_func_0_0(6) + occ_func_0_0(12) * occ_func_0_0(39) + occ_func_0_0(1) * occ_func_0_0(9) + occ_func_0_0(22) * occ_func_0_0(4) + occ_func_0_0(42) * occ_func_0_0(11) + occ_func_0_0(12) * occ_func_0_0(2) + occ_func_0_0(1) * occ_func_0_0(19) + occ_func_0_0(9) * occ_func_0_0(32) + occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(29) * occ_func_0_0(6) + occ_func_0_0(35) * occ_func_0_0(10) + occ_func_0_0(8) * occ_func_0_0(26) + occ_func_0_0(5) * occ_func_0_0(3) + occ_func_0_0(12) * occ_func_0_0(37) + occ_func_0_0(1) * occ_func_0_0(7) + occ_func_0_0(24) * occ_func_0_0(6) + occ_func_0_0(7) * occ_func_0_0(32) + occ_func_0_0(6) * occ_func_0_0(9) + occ_func_0_0(29) * occ_func_0_0(4) + occ_func_0_0(9) * occ_func_0_0(39) + occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(22) * occ_func_0_0(1) + occ_func_0_0(7) * occ_func_0_0(11) + occ_func_0_0(6) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(2) + occ_func_0_0(8) * occ_func_0_0(10) + occ_func_0_0(5) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(3) + occ_func_0_0(42) * occ_func_0_0(12) + occ_func_0_0(11) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(19) + occ_func_0_0(10) * occ_func_0_0(41) + occ_func_0_0(3) * occ_func_0_0(12) + occ_func_0_0(20) * occ_func_0_0(1) + occ_func_0_0(31) * occ_func_0_0(8) + occ_func_0_0(7) * occ_func_0_0(30) + occ_func_0_0(6) * occ_func_0_0(5) + occ_func_0_0(38) * occ_func_0_0(11) + occ_func_0_0(8) * occ_func_0_0(2) + occ_func_0_0(5) * occ_func_0_0(23) + occ_func_0_0(11) * occ_func_0_0(40) + occ_func_0_0(2) * occ_func_0_0(10) + occ_func_0_0(21) * occ_func_0_0(3) + occ_func_0_0(34) * occ_func_0_0(10) + occ_func_0_0(9) * occ_func_0_0(27) + occ_func_0_0(4) * occ_func_0_0(3) + occ_func_0_0(8) * occ_func_0_0(33) + occ_func_0_0(5) * occ_func_0_0(9) + occ_func_0_0(28) * occ_func_0_0(4)) / 24.;
  }

  /**** Basis functions for orbit 26****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 -1.0000000 1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_26_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(17) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(42) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(16) + occ_func_0_0(18) * occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(32) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(18) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(36) + occ_func_0_0(18) * occ_func_0_0(42) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(33) * occ_func_0_0(0) + occ_func_0_0(16) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(35) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(40) + occ_func_0_0(18) * occ_func_0_0(38) * occ_func_0_0(0) + occ_func_0_0(16) * occ_func_0_0(34) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(37) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(18) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(31) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(17) * occ_func_0_0(41) * occ_func_0_0(0)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_26_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(17) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(35) + occ_func_0_0(14) * occ_func_0_0(26) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(42) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(32) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(29) + occ_func_0_0(18) * occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(13) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(22) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(31) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(30) + occ_func_0_0(16) * occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(15) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(25) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(32) + occ_func_0_0(16) * occ_func_0_0(29) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(16) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(34) + occ_func_0_0(15) * occ_func_0_0(27) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(18) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(38) + occ_func_0_0(13) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(36) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(25) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(42) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(18) * occ_func_0_0(33) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(28) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(15) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(17) * occ_func_0_0(35) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(26) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(14) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(40) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(38) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(34) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(27) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(15) + occ_func_0_0(17) * occ_func_0_0(37) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(14) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(18) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(33) + occ_func_0_0(13) * occ_func_0_0(28) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(31) + occ_func_0_0(17) * occ_func_0_0(30) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(14) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(41) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(14) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(20)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_26_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(9) * occ_func_0_0(17) + occ_func_0_0(4) * occ_func_0_0(35) + occ_func_0_0(14) * occ_func_0_0(26) + occ_func_0_0(10) * occ_func_0_0(42) + occ_func_0_0(3) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(19) + occ_func_0_0(11) * occ_func_0_0(41) + occ_func_0_0(2) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(20) + occ_func_0_0(12) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(32) + occ_func_0_0(1) * occ_func_0_0(29) + occ_func_0_0(18) * occ_func_0_0(39) + occ_func_0_0(7) * occ_func_0_0(13) + occ_func_0_0(6) * occ_func_0_0(22) + occ_func_0_0(11) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(31) + occ_func_0_0(2) * occ_func_0_0(30) + occ_func_0_0(16) * occ_func_0_0(36) + occ_func_0_0(8) * occ_func_0_0(15) + occ_func_0_0(5) * occ_func_0_0(25) + occ_func_0_0(7) * occ_func_0_0(39) + occ_func_0_0(6) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(22) + occ_func_0_0(12) * occ_func_0_0(32) + occ_func_0_0(16) * occ_func_0_0(29) + occ_func_0_0(1) * occ_func_0_0(15) + occ_func_0_0(8) * occ_func_0_0(16) + occ_func_0_0(5) * occ_func_0_0(34) + occ_func_0_0(15) * occ_func_0_0(27) + occ_func_0_0(7) * occ_func_0_0(18) + occ_func_0_0(6) * occ_func_0_0(38) + occ_func_0_0(13) * occ_func_0_0(23) + occ_func_0_0(8) * occ_func_0_0(36) + occ_func_0_0(5) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(25) + occ_func_0_0(18) * occ_func_0_0(42) + occ_func_0_0(10) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(19) + occ_func_0_0(18) * occ_func_0_0(33) + occ_func_0_0(10) * occ_func_0_0(28) + occ_func_0_0(3) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(40) + occ_func_0_0(12) * occ_func_0_0(15) + occ_func_0_0(1) * occ_func_0_0(21) + occ_func_0_0(17) * occ_func_0_0(35) + occ_func_0_0(9) * occ_func_0_0(26) + occ_func_0_0(4) * occ_func_0_0(14) + occ_func_0_0(12) * occ_func_0_0(40) + occ_func_0_0(1) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(21) + occ_func_0_0(18) * occ_func_0_0(38) + occ_func_0_0(7) * occ_func_0_0(23) + occ_func_0_0(6) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(34) + occ_func_0_0(8) * occ_func_0_0(27) + occ_func_0_0(5) * occ_func_0_0(15) + occ_func_0_0(17) * occ_func_0_0(37) + occ_func_0_0(9) * occ_func_0_0(14) + occ_func_0_0(4) * occ_func_0_0(24) + occ_func_0_0(10) * occ_func_0_0(18) + occ_func_0_0(3) * occ_func_0_0(33) + occ_func_0_0(13) * occ_func_0_0(28) + occ_func_0_0(11) * occ_func_0_0(31) + occ_func_0_0(17) * occ_func_0_0(30) + occ_func_0_0(2) * occ_func_0_0(14) + occ_func_0_0(9) * occ_func_0_0(37) + occ_func_0_0(4) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(24) + occ_func_0_0(17) * occ_func_0_0(41) + occ_func_0_0(11) * occ_func_0_0(14) + occ_func_0_0(2) * occ_func_0_0(20)) / 24.;
  }

  /**** Basis functions for orbit 27****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 1.0000000 A  B  

1.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_27_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(32) * occ_func_0_0(39) + occ_func_0_0(35) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(38) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(0) * occ_func_0_0(41) * occ_func_0_0(42) + occ_func_0_0(36) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(40) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(36) + occ_func_0_0(42) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(0) * occ_func_0_0(38) * occ_func_0_0(31) + occ_func_0_0(0) * occ_func_0_0(37) * occ_func_0_0(32) + occ_func_0_0(33) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(35) * occ_func_0_0(34) + occ_func_0_0(41) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(0) * occ_func_0_0(39) * occ_func_0_0(37) + occ_func_0_0(39) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(39) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(37) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(38) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(34) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(40) * occ_func_0_0(41) + occ_func_0_0(42) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(33) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(36) * occ_func_0_0(38)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_27_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(32) * occ_func_0_0(39) + occ_func_0_0(29) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(22) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(35) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(34) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(26) + occ_func_0_0(38) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(36) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(23) * occ_func_0_0(25) + occ_func_0_0(0) * occ_func_0_0(41) * occ_func_0_0(42) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(19) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(36) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(31) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25) * occ_func_0_0(30) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(40) + occ_func_0_0(21) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(33) + occ_func_0_0(28) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(36) + occ_func_0_0(30) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(25) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(41) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(20) + occ_func_0_0(0) * occ_func_0_0(38) * occ_func_0_0(31) + occ_func_0_0(30) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(23) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(37) * occ_func_0_0(32) + occ_func_0_0(29) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(24) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(33) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(35) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(26) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(41) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(40) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(21) + occ_func_0_0(0) * occ_func_0_0(39) * occ_func_0_0(37) + occ_func_0_0(24) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(22) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(39) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(32) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(22) * occ_func_0_0(29) + occ_func_0_0(39) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(37) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(24) * occ_func_0_0(22) + occ_func_0_0(37) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(32) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(29) * occ_func_0_0(24) + occ_func_0_0(38) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(31) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(30) * occ_func_0_0(23) + occ_func_0_0(34) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(33) * occ_func_0_0(27) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(40) * occ_func_0_0(41) + occ_func_0_0(21) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(20) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(40) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(0) * occ_func_0_0(21) * occ_func_0_0(19) + occ_func_0_0(33) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(34) * occ_func_0_0(28) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(36) * occ_func_0_0(38) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(23) * occ_func_0_0(4) * occ_func_0_0(0)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_27_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(32) * occ_func_0_0(39) + occ_func_0_0(29) * occ_func_0_0(11) + occ_func_0_0(22) * occ_func_0_0(2) + occ_func_0_0(35) * occ_func_0_0(7) + occ_func_0_0(34) * occ_func_0_0(6) + occ_func_0_0(27) * occ_func_0_0(26) + occ_func_0_0(38) * occ_func_0_0(9) + occ_func_0_0(36) * occ_func_0_0(4) + occ_func_0_0(23) * occ_func_0_0(25) + occ_func_0_0(41) * occ_func_0_0(42) + occ_func_0_0(20) * occ_func_0_0(8) + occ_func_0_0(19) * occ_func_0_0(5) + occ_func_0_0(36) * occ_func_0_0(10) + occ_func_0_0(31) * occ_func_0_0(3) + occ_func_0_0(25) * occ_func_0_0(30) + occ_func_0_0(42) * occ_func_0_0(40) + occ_func_0_0(21) * occ_func_0_0(9) + occ_func_0_0(19) * occ_func_0_0(4) + occ_func_0_0(35) * occ_func_0_0(12) + occ_func_0_0(26) * occ_func_0_0(33) + occ_func_0_0(28) * occ_func_0_0(1) + occ_func_0_0(31) * occ_func_0_0(36) + occ_func_0_0(30) * occ_func_0_0(10) + occ_func_0_0(25) * occ_func_0_0(3) + occ_func_0_0(42) * occ_func_0_0(8) + occ_func_0_0(41) * occ_func_0_0(5) + occ_func_0_0(19) * occ_func_0_0(20) + occ_func_0_0(38) * occ_func_0_0(31) + occ_func_0_0(30) * occ_func_0_0(12) + occ_func_0_0(23) * occ_func_0_0(1) + occ_func_0_0(37) * occ_func_0_0(32) + occ_func_0_0(29) * occ_func_0_0(10) + occ_func_0_0(24) * occ_func_0_0(3) + occ_func_0_0(33) * occ_func_0_0(12) + occ_func_0_0(28) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(1) + occ_func_0_0(35) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(7) + occ_func_0_0(26) * occ_func_0_0(6) + occ_func_0_0(41) * occ_func_0_0(7) + occ_func_0_0(40) * occ_func_0_0(6) + occ_func_0_0(20) * occ_func_0_0(21) + occ_func_0_0(39) * occ_func_0_0(37) + occ_func_0_0(24) * occ_func_0_0(8) + occ_func_0_0(22) * occ_func_0_0(5) + occ_func_0_0(39) * occ_func_0_0(11) + occ_func_0_0(32) * occ_func_0_0(2) + occ_func_0_0(22) * occ_func_0_0(29) + occ_func_0_0(39) * occ_func_0_0(8) + occ_func_0_0(37) * occ_func_0_0(5) + occ_func_0_0(24) * occ_func_0_0(22) + occ_func_0_0(37) * occ_func_0_0(10) + occ_func_0_0(32) * occ_func_0_0(3) + occ_func_0_0(29) * occ_func_0_0(24) + occ_func_0_0(38) * occ_func_0_0(12) + occ_func_0_0(31) * occ_func_0_0(1) + occ_func_0_0(30) * occ_func_0_0(23) + occ_func_0_0(34) * occ_func_0_0(11) + occ_func_0_0(33) * occ_func_0_0(27) + occ_func_0_0(28) * occ_func_0_0(2) + occ_func_0_0(40) * occ_func_0_0(41) + occ_func_0_0(21) * occ_func_0_0(7) + occ_func_0_0(20) * occ_func_0_0(6) + occ_func_0_0(42) * occ_func_0_0(9) + occ_func_0_0(40) * occ_func_0_0(4) + occ_func_0_0(21) * occ_func_0_0(19) + occ_func_0_0(33) * occ_func_0_0(11) + occ_func_0_0(34) * occ_func_0_0(28) + occ_func_0_0(27) * occ_func_0_0(2) + occ_func_0_0(36) * occ_func_0_0(38) + occ_func_0_0(25) * occ_func_0_0(9) + occ_func_0_0(23) * occ_func_0_0(4)) / 24.;
  }

  /**** Basis functions for orbit 28****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 -2.0000000 A  B  

1.0000000 -1.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_28_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(16) + occ_func_0_0(39) * occ_func_0_0(0) * occ_func_0_0(42) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(17) * occ_func_0_0(32) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(16) + occ_func_0_0(18) * occ_func_0_0(31) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(40) * occ_func_0_0(17) + occ_func_0_0(33) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(0) * occ_func_0_0(36) * occ_func_0_0(18) + occ_func_0_0(40) * occ_func_0_0(0) * occ_func_0_0(34) + occ_func_0_0(42) * occ_func_0_0(0) * occ_func_0_0(38) + occ_func_0_0(0) * occ_func_0_0(32) * occ_func_0_0(18) + occ_func_0_0(17) * occ_func_0_0(42) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(39) * occ_func_0_0(17) + occ_func_0_0(18) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(41) * occ_func_0_0(18) + occ_func_0_0(36) * occ_func_0_0(0) * occ_func_0_0(40) + occ_func_0_0(16) * occ_func_0_0(38) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(34) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(37) * occ_func_0_0(0) + occ_func_0_0(38) * occ_func_0_0(0) * occ_func_0_0(33) + occ_func_0_0(0) * occ_func_0_0(35) * occ_func_0_0(16) + occ_func_0_0(41) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(16) * occ_func_0_0(41) * occ_func_0_0(0)) / 24.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_28_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(16) + occ_func_0_0(30) * occ_func_0_0(0) * occ_func_0_0(35) + occ_func_0_0(15) * occ_func_0_0(26) * occ_func_0_0(0) + occ_func_0_0(39) * occ_func_0_0(0) * occ_func_0_0(42) + occ_func_0_0(0) * occ_func_0_0(22) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(32) * occ_func_0_0(0) + occ_func_0_0(34) * occ_func_0_0(0) * occ_func_0_0(29) + occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(14) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(22) + occ_func_0_0(18) * occ_func_0_0(31) * occ_func_0_0(0) + occ_func_0_0(37) * occ_func_0_0(0) * occ_func_0_0(30) + occ_func_0_0(0) * occ_func_0_0(24) * occ_func_0_0(13) + occ_func_0_0(0) * occ_func_0_0(40) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(21) * occ_func_0_0(0) * occ_func_0_0(25) + occ_func_0_0(33) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(36) * occ_func_0_0(18) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(32) + occ_func_0_0(13) * occ_func_0_0(29) * occ_func_0_0(0) + occ_func_0_0(40) * occ_func_0_0(0) * occ_func_0_0(34) + occ_func_0_0(18) * occ_func_0_0(27) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(21) * occ_func_0_0(13) + occ_func_0_0(42) * occ_func_0_0(0) * occ_func_0_0(38) + occ_func_0_0(17) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(14) + occ_func_0_0(0) * occ_func_0_0(32) * occ_func_0_0(18) + occ_func_0_0(29) * occ_func_0_0(0) * occ_func_0_0(36) + occ_func_0_0(13) * occ_func_0_0(25) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(42) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(38) * occ_func_0_0(14) + occ_func_0_0(23) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(39) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(33) * occ_func_0_0(0) + occ_func_0_0(22) * occ_func_0_0(0) * occ_func_0_0(28) + occ_func_0_0(18) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(34) * occ_func_0_0(13) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(0) * occ_func_0_0(41) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(35) * occ_func_0_0(0) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(26) + occ_func_0_0(36) * occ_func_0_0(0) * occ_func_0_0(40) + occ_func_0_0(0) * occ_func_0_0(25) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(16) * occ_func_0_0(38) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(33) * occ_func_0_0(15) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(17) * occ_func_0_0(34) * occ_func_0_0(0) + occ_func_0_0(32) * occ_func_0_0(0) * occ_func_0_0(27) + occ_func_0_0(0) * occ_func_0_0(29) * occ_func_0_0(14) + occ_func_0_0(18) * occ_func_0_0(37) * occ_func_0_0(0) + occ_func_0_0(31) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(0) * occ_func_0_0(30) * occ_func_0_0(13) + occ_func_0_0(38) * occ_func_0_0(0) * occ_func_0_0(33) + occ_func_0_0(16) * occ_func_0_0(28) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(23) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(35) * occ_func_0_0(16) + occ_func_0_0(26) * occ_func_0_0(0) * occ_func_0_0(31) + occ_func_0_0(15) * occ_func_0_0(30) * occ_func_0_0(0) + occ_func_0_0(41) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(16) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(15) + occ_func_0_0(16) * occ_func_0_0(41) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(37) * occ_func_0_0(15) + occ_func_0_0(24) * occ_func_0_0(0) * occ_func_0_0(20)) / 24.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_28_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(31) * occ_func_0_0(16) + occ_func_0_0(30) * occ_func_0_0(35) + occ_func_0_0(15) * occ_func_0_0(26) + occ_func_0_0(39) * occ_func_0_0(42) + occ_func_0_0(22) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(19) + occ_func_0_0(35) * occ_func_0_0(41) + occ_func_0_0(26) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(20) + occ_func_0_0(17) * occ_func_0_0(32) + occ_func_0_0(34) * occ_func_0_0(29) + occ_func_0_0(27) * occ_func_0_0(14) + occ_func_0_0(42) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(39) + occ_func_0_0(19) * occ_func_0_0(22) + occ_func_0_0(18) * occ_func_0_0(31) + occ_func_0_0(37) * occ_func_0_0(30) + occ_func_0_0(24) * occ_func_0_0(13) + occ_func_0_0(40) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(36) + occ_func_0_0(21) * occ_func_0_0(25) + occ_func_0_0(33) * occ_func_0_0(39) + occ_func_0_0(28) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(22) + occ_func_0_0(36) * occ_func_0_0(18) + occ_func_0_0(25) * occ_func_0_0(32) + occ_func_0_0(13) * occ_func_0_0(29) + occ_func_0_0(40) * occ_func_0_0(34) + occ_func_0_0(18) * occ_func_0_0(27) + occ_func_0_0(21) * occ_func_0_0(13) + occ_func_0_0(42) * occ_func_0_0(38) + occ_func_0_0(17) * occ_func_0_0(23) + occ_func_0_0(19) * occ_func_0_0(14) + occ_func_0_0(32) * occ_func_0_0(18) + occ_func_0_0(29) * occ_func_0_0(36) + occ_func_0_0(13) * occ_func_0_0(25) + occ_func_0_0(17) * occ_func_0_0(42) + occ_func_0_0(38) * occ_func_0_0(14) + occ_func_0_0(23) * occ_func_0_0(19) + occ_func_0_0(39) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(33) + occ_func_0_0(22) * occ_func_0_0(28) + occ_func_0_0(18) * occ_func_0_0(40) + occ_func_0_0(34) * occ_func_0_0(13) + occ_func_0_0(27) * occ_func_0_0(21) + occ_func_0_0(41) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(35) + occ_func_0_0(20) * occ_func_0_0(26) + occ_func_0_0(36) * occ_func_0_0(40) + occ_func_0_0(25) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(21) + occ_func_0_0(16) * occ_func_0_0(38) + occ_func_0_0(33) * occ_func_0_0(15) + occ_func_0_0(28) * occ_func_0_0(23) + occ_func_0_0(17) * occ_func_0_0(34) + occ_func_0_0(32) * occ_func_0_0(27) + occ_func_0_0(29) * occ_func_0_0(14) + occ_func_0_0(18) * occ_func_0_0(37) + occ_func_0_0(31) * occ_func_0_0(24) + occ_func_0_0(30) * occ_func_0_0(13) + occ_func_0_0(38) * occ_func_0_0(33) + occ_func_0_0(16) * occ_func_0_0(28) + occ_func_0_0(23) * occ_func_0_0(15) + occ_func_0_0(35) * occ_func_0_0(16) + occ_func_0_0(26) * occ_func_0_0(31) + occ_func_0_0(15) * occ_func_0_0(30) + occ_func_0_0(41) * occ_func_0_0(37) + occ_func_0_0(16) * occ_func_0_0(24) + occ_func_0_0(20) * occ_func_0_0(15) + occ_func_0_0(16) * occ_func_0_0(41) + occ_func_0_0(37) * occ_func_0_0(15) + occ_func_0_0(24) * occ_func_0_0(20)) / 24.;
  }

  /**** Basis functions for orbit 29****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 -2.0000000 A  B  

0.0000000 2.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_29_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(33) + occ_func_0_0(41) * occ_func_0_0(34) * occ_func_0_0(0) + occ_func_0_0(40) * occ_func_0_0(38) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(37) + occ_func_0_0(41) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(0) * occ_func_0_0(36) * occ_func_0_0(42) + occ_func_0_0(33) * occ_func_0_0(0) * occ_func_0_0(32) + occ_func_0_0(35) * occ_func_0_0(40) * occ_func_0_0(0)) / 8.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_29_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(33) + occ_func_0_0(30) * occ_func_0_0(0) * occ_func_0_0(32) + occ_func_0_0(28) * occ_func_0_0(29) * occ_func_0_0(0) + occ_func_0_0(41) * occ_func_0_0(34) * occ_func_0_0(0) + occ_func_0_0(39) * occ_func_0_0(0) * occ_func_0_0(27) + occ_func_0_0(0) * occ_func_0_0(22) * occ_func_0_0(20) + occ_func_0_0(40) * occ_func_0_0(38) * occ_func_0_0(0) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(21) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(37) + occ_func_0_0(24) * occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(25) + occ_func_0_0(41) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(34) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(27) + occ_func_0_0(0) * occ_func_0_0(36) * occ_func_0_0(42) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(19) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(33) * occ_func_0_0(0) * occ_func_0_0(32) + occ_func_0_0(31) * occ_func_0_0(29) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(30) + occ_func_0_0(35) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(38) * occ_func_0_0(26) + occ_func_0_0(23) * occ_func_0_0(0) * occ_func_0_0(21)) / 8.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_29_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(31) * occ_func_0_0(33) + occ_func_0_0(30) * occ_func_0_0(32) + occ_func_0_0(28) * occ_func_0_0(29) + occ_func_0_0(41) * occ_func_0_0(34) + occ_func_0_0(39) * occ_func_0_0(27) + occ_func_0_0(22) * occ_func_0_0(20) + occ_func_0_0(40) * occ_func_0_0(38) + occ_func_0_0(35) * occ_func_0_0(23) + occ_func_0_0(26) * occ_func_0_0(21) + occ_func_0_0(42) * occ_func_0_0(37) + occ_func_0_0(24) * occ_func_0_0(36) + occ_func_0_0(19) * occ_func_0_0(25) + occ_func_0_0(41) * occ_func_0_0(39) + occ_func_0_0(34) * occ_func_0_0(22) + occ_func_0_0(20) * occ_func_0_0(27) + occ_func_0_0(36) * occ_func_0_0(42) + occ_func_0_0(25) * occ_func_0_0(37) + occ_func_0_0(19) * occ_func_0_0(24) + occ_func_0_0(33) * occ_func_0_0(32) + occ_func_0_0(31) * occ_func_0_0(29) + occ_func_0_0(28) * occ_func_0_0(30) + occ_func_0_0(35) * occ_func_0_0(40) + occ_func_0_0(38) * occ_func_0_0(26) + occ_func_0_0(23) * occ_func_0_0(21)) / 8.;
  }

  /**** Basis functions for orbit 30****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_30_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(51) + occ_func_0_0(52) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(53) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(49) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(50)) / 6.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_30_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(51) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(46) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(52) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(45) + occ_func_0_0(53) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(44) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(54) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(43) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(49) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(48) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(50) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(47) * occ_func_0_0(5) * occ_func_0_0(0)) / 6.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_30_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(9) * occ_func_0_0(51) + occ_func_0_0(4) * occ_func_0_0(9) + occ_func_0_0(46) * occ_func_0_0(4) + occ_func_0_0(52) * occ_func_0_0(10) + occ_func_0_0(10) * occ_func_0_0(3) + occ_func_0_0(3) * occ_func_0_0(45) + occ_func_0_0(53) * occ_func_0_0(11) + occ_func_0_0(11) * occ_func_0_0(2) + occ_func_0_0(2) * occ_func_0_0(44) + occ_func_0_0(12) * occ_func_0_0(54) + occ_func_0_0(1) * occ_func_0_0(12) + occ_func_0_0(43) * occ_func_0_0(1) + occ_func_0_0(7) * occ_func_0_0(49) + occ_func_0_0(6) * occ_func_0_0(7) + occ_func_0_0(48) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(50) + occ_func_0_0(5) * occ_func_0_0(8) + occ_func_0_0(47) * occ_func_0_0(5)) / 6.;
  }

  /**** Basis functions for orbit 31****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

1.0000000 1.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_31_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_0(18) + occ_func_0_0(52) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_0(18) + occ_func_0_0(17) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(16) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(52) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(18) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_31_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_0(18) + occ_func_0_0(46) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(13) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(52) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(17) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(45) * occ_func_0_0(14) + occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(16) * occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(43) * occ_func_0_0(0) * occ_func_0_0(14) + occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(17) + occ_func_0_0(48) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(14) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(44) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_0(18) + occ_func_0_0(47) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(13) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(50) * occ_func_0_0(14) + occ_func_0_0(47) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(15) + occ_func_0_0(48) * occ_func_0_0(0) * occ_func_0_0(14) + occ_func_0_0(52) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(45) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(54) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(43) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_0(15) + occ_func_0_0(46) * occ_func_0_0(0) * occ_func_0_0(13)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_31_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(51) * occ_func_0_0(18) + occ_func_0_0(46) * occ_func_0_0(16) + occ_func_0_0(13) * occ_func_0_0(15) + occ_func_0_0(52) * occ_func_0_0(16) + occ_func_0_0(17) * occ_func_0_0(15) + occ_func_0_0(45) * occ_func_0_0(14) + occ_func_0_0(53) * occ_func_0_0(18) + occ_func_0_0(16) * occ_func_0_0(13) + occ_func_0_0(44) * occ_func_0_0(15) + occ_func_0_0(54) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(17) + occ_func_0_0(43) * occ_func_0_0(14) + occ_func_0_0(49) * occ_func_0_0(17) + occ_func_0_0(48) * occ_func_0_0(16) + occ_func_0_0(14) * occ_func_0_0(15) + occ_func_0_0(53) * occ_func_0_0(16) + occ_func_0_0(15) * occ_func_0_0(18) + occ_func_0_0(44) * occ_func_0_0(13) + occ_func_0_0(50) * occ_func_0_0(18) + occ_func_0_0(47) * occ_func_0_0(17) + occ_func_0_0(13) * occ_func_0_0(14) + occ_func_0_0(17) * occ_func_0_0(18) + occ_func_0_0(50) * occ_func_0_0(14) + occ_func_0_0(47) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(17) + occ_func_0_0(49) * occ_func_0_0(15) + occ_func_0_0(48) * occ_func_0_0(14) + occ_func_0_0(52) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(14) + occ_func_0_0(45) * occ_func_0_0(15) + occ_func_0_0(54) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(18) + occ_func_0_0(43) * occ_func_0_0(13) + occ_func_0_0(16) * occ_func_0_0(18) + occ_func_0_0(51) * occ_func_0_0(15) + occ_func_0_0(46) * occ_func_0_0(13)) / 12.;
  }

  /**** Basis functions for orbit 32****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 0.0000000 1.0000000 A  B  

0.0000000 2.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_32_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(33) + occ_func_0_0(41) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(40) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(39) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(36) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(32) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(37) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(35) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(36) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(41) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(52) + occ_func_0_0(42) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(33) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(32) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(35) + occ_func_0_0(38) * occ_func_0_0(50) * occ_func_0_0(0) + occ_func_0_0(39) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(34) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(42) + occ_func_0_0(31) * occ_func_0_0(50) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(32) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(38) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(42) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(37) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(52) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(31) + occ_func_0_0(40) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(34) + occ_func_0_0(41) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(40) + occ_func_0_0(39) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(37) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(38) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(33) * occ_func_0_0(50) * occ_func_0_0(0)) / 48.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_32_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(33) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(28) * occ_func_0_0(47) * occ_func_0_0(0) + occ_func_0_0(41) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(43) + occ_func_0_0(40) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(21) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(39) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(22) * occ_func_0_0(46) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(37) + occ_func_0_0(24) * occ_func_0_0(43) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(48) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(38) + occ_func_0_0(23) * occ_func_0_0(44) * occ_func_0_0(0) + occ_func_0_0(32) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(29) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(46) + occ_func_0_0(37) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(48) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(24) + occ_func_0_0(35) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(26) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(36) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(25) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(44) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(31) * occ_func_0_0(48) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(30) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(40) + occ_func_0_0(21) * occ_func_0_0(44) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(34) * occ_func_0_0(46) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(27) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(41) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(52) + occ_func_0_0(20) * occ_func_0_0(45) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(52) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(45) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(44) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(19) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(31) + occ_func_0_0(30) * occ_func_0_0(47) * occ_func_0_0(0) + occ_func_0_0(33) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(28) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(46) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(32) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(29) * occ_func_0_0(48) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(50) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(47) + occ_func_0_0(38) * occ_func_0_0(50) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(47) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(23) + occ_func_0_0(39) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(22) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(43) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(42) + occ_func_0_0(19) * occ_func_0_0(43) * occ_func_0_0(0) + occ_func_0_0(34) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(27) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(42) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(19) * occ_func_0_0(44) * occ_func_0_0(0) + occ_func_0_0(31) * occ_func_0_0(50) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(30) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(47) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(33) + occ_func_0_0(28) * occ_func_0_0(46) * occ_func_0_0(0) + occ_func_0_0(32) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(48) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(29) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(35) * occ_func_0_0(47) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(26) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(38) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(23) * occ_func_0_0(47) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(39) + occ_func_0_0(22) * occ_func_0_0(43) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(43) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(32) + occ_func_0_0(29) * occ_func_0_0(46) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(37) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(24) * occ_func_0_0(48) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(52) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(45) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(53) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(44) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(31) + occ_func_0_0(30) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(48) + occ_func_0_0(40) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(21) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(44) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(46) + occ_func_0_0(41) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(20) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(41) + occ_func_0_0(20) * occ_func_0_0(43) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(40) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(52) + occ_func_0_0(21) * occ_func_0_0(45) * occ_func_0_0(0) + occ_func_0_0(39) * occ_func_0_0(51) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(46) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(22) + occ_func_0_0(37) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(24) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(43) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(49) + occ_func_0_0(36) * occ_func_0_0(48) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(25) + occ_func_0_0(38) * occ_func_0_0(53) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(23) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(44) + occ_func_0_0(33) * occ_func_0_0(50) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(47) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(28)) / 48.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_32_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(33) + occ_func_0_0(6) * occ_func_0_0(50) + occ_func_0_0(28) * occ_func_0_0(47) + occ_func_0_0(41) * occ_func_0_0(54) + occ_func_0_0(9) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(43) + occ_func_0_0(40) * occ_func_0_0(52) + occ_func_0_0(8) * occ_func_0_0(45) + occ_func_0_0(5) * occ_func_0_0(21) + occ_func_0_0(10) * occ_func_0_0(39) + occ_func_0_0(3) * occ_func_0_0(51) + occ_func_0_0(22) * occ_func_0_0(46) + occ_func_0_0(11) * occ_func_0_0(54) + occ_func_0_0(2) * occ_func_0_0(37) + occ_func_0_0(24) * occ_func_0_0(43) + occ_func_0_0(12) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(49) + occ_func_0_0(1) * occ_func_0_0(48) + occ_func_0_0(10) * occ_func_0_0(53) + occ_func_0_0(3) * occ_func_0_0(38) + occ_func_0_0(23) * occ_func_0_0(44) + occ_func_0_0(32) * occ_func_0_0(51) + occ_func_0_0(8) * occ_func_0_0(29) + occ_func_0_0(5) * occ_func_0_0(46) + occ_func_0_0(37) * occ_func_0_0(49) + occ_func_0_0(11) * occ_func_0_0(48) + occ_func_0_0(2) * occ_func_0_0(24) + occ_func_0_0(35) * occ_func_0_0(52) + occ_func_0_0(11) * occ_func_0_0(26) + occ_func_0_0(2) * occ_func_0_0(45) + occ_func_0_0(36) * occ_func_0_0(53) + occ_func_0_0(12) * occ_func_0_0(25) + occ_func_0_0(1) * occ_func_0_0(44) + occ_func_0_0(9) * occ_func_0_0(49) + occ_func_0_0(31) * occ_func_0_0(48) + occ_func_0_0(4) * occ_func_0_0(30) + occ_func_0_0(8) * occ_func_0_0(53) + occ_func_0_0(5) * occ_func_0_0(40) + occ_func_0_0(21) * occ_func_0_0(44) + occ_func_0_0(12) * occ_func_0_0(51) + occ_func_0_0(34) * occ_func_0_0(46) + occ_func_0_0(1) * occ_func_0_0(27) + occ_func_0_0(9) * occ_func_0_0(41) + occ_func_0_0(4) * occ_func_0_0(52) + occ_func_0_0(20) * occ_func_0_0(45) + occ_func_0_0(12) * occ_func_0_0(52) + occ_func_0_0(1) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(45) + occ_func_0_0(42) * occ_func_0_0(53) + occ_func_0_0(7) * occ_func_0_0(44) + occ_func_0_0(6) * occ_func_0_0(19) + occ_func_0_0(9) * occ_func_0_0(50) + occ_func_0_0(4) * occ_func_0_0(31) + occ_func_0_0(30) * occ_func_0_0(47) + occ_func_0_0(33) * occ_func_0_0(51) + occ_func_0_0(7) * occ_func_0_0(28) + occ_func_0_0(6) * occ_func_0_0(46) + occ_func_0_0(8) * occ_func_0_0(32) + occ_func_0_0(5) * occ_func_0_0(49) + occ_func_0_0(29) * occ_func_0_0(48) + occ_func_0_0(11) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(50) + occ_func_0_0(2) * occ_func_0_0(47) + occ_func_0_0(38) * occ_func_0_0(50) + occ_func_0_0(10) * occ_func_0_0(47) + occ_func_0_0(3) * occ_func_0_0(23) + occ_func_0_0(39) * occ_func_0_0(54) + occ_func_0_0(10) * occ_func_0_0(22) + occ_func_0_0(3) * occ_func_0_0(43) + occ_func_0_0(7) * occ_func_0_0(54) + occ_func_0_0(6) * occ_func_0_0(42) + occ_func_0_0(19) * occ_func_0_0(43) + occ_func_0_0(34) * occ_func_0_0(52) + occ_func_0_0(12) * occ_func_0_0(27) + occ_func_0_0(1) * occ_func_0_0(45) + occ_func_0_0(7) * occ_func_0_0(42) + occ_func_0_0(6) * occ_func_0_0(53) + occ_func_0_0(19) * occ_func_0_0(44) + occ_func_0_0(31) * occ_func_0_0(50) + occ_func_0_0(9) * occ_func_0_0(30) + occ_func_0_0(4) * occ_func_0_0(47) + occ_func_0_0(7) * occ_func_0_0(51) + occ_func_0_0(6) * occ_func_0_0(33) + occ_func_0_0(28) * occ_func_0_0(46) + occ_func_0_0(32) * occ_func_0_0(49) + occ_func_0_0(8) * occ_func_0_0(48) + occ_func_0_0(5) * occ_func_0_0(29) + occ_func_0_0(11) * occ_func_0_0(50) + occ_func_0_0(35) * occ_func_0_0(47) + occ_func_0_0(2) * occ_func_0_0(26) + occ_func_0_0(10) * occ_func_0_0(38) + occ_func_0_0(3) * occ_func_0_0(50) + occ_func_0_0(23) * occ_func_0_0(47) + occ_func_0_0(10) * occ_func_0_0(54) + occ_func_0_0(3) * occ_func_0_0(39) + occ_func_0_0(22) * occ_func_0_0(43) + occ_func_0_0(42) * occ_func_0_0(54) + occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(6) * occ_func_0_0(43) + occ_func_0_0(8) * occ_func_0_0(51) + occ_func_0_0(5) * occ_func_0_0(32) + occ_func_0_0(29) * occ_func_0_0(46) + occ_func_0_0(11) * occ_func_0_0(37) + occ_func_0_0(2) * occ_func_0_0(49) + occ_func_0_0(24) * occ_func_0_0(48) + occ_func_0_0(11) * occ_func_0_0(52) + occ_func_0_0(2) * occ_func_0_0(35) + occ_func_0_0(26) * occ_func_0_0(45) + occ_func_0_0(12) * occ_func_0_0(53) + occ_func_0_0(1) * occ_func_0_0(36) + occ_func_0_0(25) * occ_func_0_0(44) + occ_func_0_0(9) * occ_func_0_0(31) + occ_func_0_0(30) * occ_func_0_0(49) + occ_func_0_0(4) * occ_func_0_0(48) + occ_func_0_0(40) * occ_func_0_0(53) + occ_func_0_0(8) * occ_func_0_0(21) + occ_func_0_0(5) * occ_func_0_0(44) + occ_func_0_0(12) * occ_func_0_0(34) + occ_func_0_0(27) * occ_func_0_0(51) + occ_func_0_0(1) * occ_func_0_0(46) + occ_func_0_0(41) * occ_func_0_0(52) + occ_func_0_0(9) * occ_func_0_0(45) + occ_func_0_0(4) * occ_func_0_0(20) + occ_func_0_0(9) * occ_func_0_0(54) + occ_func_0_0(4) * occ_func_0_0(41) + occ_func_0_0(20) * occ_func_0_0(43) + occ_func_0_0(8) * occ_func_0_0(40) + occ_func_0_0(5) * occ_func_0_0(52) + occ_func_0_0(21) * occ_func_0_0(45) + occ_func_0_0(39) * occ_func_0_0(51) + occ_func_0_0(10) * occ_func_0_0(46) + occ_func_0_0(3) * occ_func_0_0(22) + occ_func_0_0(37) * occ_func_0_0(54) + occ_func_0_0(11) * occ_func_0_0(24) + occ_func_0_0(2) * occ_func_0_0(43) + occ_func_0_0(12) * occ_func_0_0(49) + occ_func_0_0(36) * occ_func_0_0(48) + occ_func_0_0(1) * occ_func_0_0(25) + occ_func_0_0(38) * occ_func_0_0(53) + occ_func_0_0(10) * occ_func_0_0(23) + occ_func_0_0(3) * occ_func_0_0(44) + occ_func_0_0(33) * occ_func_0_0(50) + occ_func_0_0(7) * occ_func_0_0(47) + occ_func_0_0(6) * occ_func_0_0(28)) / 48.;
  }

  /**** Basis functions for orbit 33****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 1.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_33_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(32) * occ_func_0_0(54) + occ_func_0_0(34) * occ_func_0_0(0) * occ_func_0_0(36) + occ_func_0_0(38) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(0) * occ_func_0_0(41) * occ_func_0_0(53) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(52) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(53) + occ_func_0_0(42) * occ_func_0_0(0) * occ_func_0_0(33) + occ_func_0_0(49) * occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(37) * occ_func_0_0(51) + occ_func_0_0(54) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(50) * occ_func_0_0(36) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_33_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(32) * occ_func_0_0(54) + occ_func_0_0(29) * occ_func_0_0(0) * occ_func_0_0(40) + occ_func_0_0(43) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(34) * occ_func_0_0(0) * occ_func_0_0(36) + occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(50) + occ_func_0_0(47) * occ_func_0_0(25) * occ_func_0_0(0) + occ_func_0_0(38) * occ_func_0_0(0) * occ_func_0_0(39) + occ_func_0_0(0) * occ_func_0_0(23) * occ_func_0_0(49) + occ_func_0_0(48) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(41) * occ_func_0_0(53) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(31) + occ_func_0_0(44) * occ_func_0_0(30) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(42) * occ_func_0_0(52) + occ_func_0_0(45) * occ_func_0_0(33) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(28) + occ_func_0_0(35) * occ_func_0_0(0) * occ_func_0_0(37) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(51) + occ_func_0_0(46) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(31) * occ_func_0_0(53) + occ_func_0_0(30) * occ_func_0_0(0) * occ_func_0_0(41) + occ_func_0_0(44) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(42) * occ_func_0_0(0) * occ_func_0_0(33) + occ_func_0_0(52) * occ_func_0_0(28) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(45) + occ_func_0_0(49) * occ_func_0_0(39) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(38) * occ_func_0_0(48) + occ_func_0_0(23) * occ_func_0_0(0) * occ_func_0_0(22) + occ_func_0_0(0) * occ_func_0_0(37) * occ_func_0_0(51) + occ_func_0_0(46) * occ_func_0_0(35) * occ_func_0_0(0) + occ_func_0_0(24) * occ_func_0_0(0) * occ_func_0_0(26) + occ_func_0_0(54) * occ_func_0_0(40) * occ_func_0_0(0) + occ_func_0_0(32) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(0) * occ_func_0_0(29) * occ_func_0_0(43) + occ_func_0_0(50) * occ_func_0_0(36) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(34) * occ_func_0_0(47) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(25)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_33_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(32) * occ_func_0_0(54) + occ_func_0_0(29) * occ_func_0_0(40) + occ_func_0_0(43) * occ_func_0_0(21) + occ_func_0_0(34) * occ_func_0_0(36) + occ_func_0_0(27) * occ_func_0_0(50) + occ_func_0_0(47) * occ_func_0_0(25) + occ_func_0_0(38) * occ_func_0_0(39) + occ_func_0_0(23) * occ_func_0_0(49) + occ_func_0_0(48) * occ_func_0_0(22) + occ_func_0_0(41) * occ_func_0_0(53) + occ_func_0_0(20) * occ_func_0_0(31) + occ_func_0_0(44) * occ_func_0_0(30) + occ_func_0_0(42) * occ_func_0_0(52) + occ_func_0_0(45) * occ_func_0_0(33) + occ_func_0_0(19) * occ_func_0_0(28) + occ_func_0_0(35) * occ_func_0_0(37) + occ_func_0_0(26) * occ_func_0_0(51) + occ_func_0_0(46) * occ_func_0_0(24) + occ_func_0_0(31) * occ_func_0_0(53) + occ_func_0_0(30) * occ_func_0_0(41) + occ_func_0_0(44) * occ_func_0_0(20) + occ_func_0_0(42) * occ_func_0_0(33) + occ_func_0_0(52) * occ_func_0_0(28) + occ_func_0_0(19) * occ_func_0_0(45) + occ_func_0_0(49) * occ_func_0_0(39) + occ_func_0_0(38) * occ_func_0_0(48) + occ_func_0_0(23) * occ_func_0_0(22) + occ_func_0_0(37) * occ_func_0_0(51) + occ_func_0_0(46) * occ_func_0_0(35) + occ_func_0_0(24) * occ_func_0_0(26) + occ_func_0_0(54) * occ_func_0_0(40) + occ_func_0_0(32) * occ_func_0_0(21) + occ_func_0_0(29) * occ_func_0_0(43) + occ_func_0_0(50) * occ_func_0_0(36) + occ_func_0_0(34) * occ_func_0_0(47) + occ_func_0_0(27) * occ_func_0_0(25)) / 12.;
  }

  /**** Basis functions for orbit 34****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 0.0000000 2.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_34_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(51) + occ_func_0_0(52) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(53) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(54) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_0(50) + occ_func_0_0(54) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(51) * occ_func_0_0(0) * occ_func_0_0(50)) / 8.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_34_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(49) * occ_func_0_0(51) + occ_func_0_0(48) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(46) * occ_func_0_0(47) * occ_func_0_0(0) + occ_func_0_0(52) * occ_func_0_0(54) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(51) * occ_func_0_0(45) + occ_func_0_0(46) * occ_func_0_0(0) * occ_func_0_0(43) + occ_func_0_0(53) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(50) * occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(0) * occ_func_0_0(47) * occ_func_0_0(44) + occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_0(54) + occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_0(49) + occ_func_0_0(48) * occ_func_0_0(43) * occ_func_0_0(0) + occ_func_0_0(54) * occ_func_0_0(49) * occ_func_0_0(0) + occ_func_0_0(53) * occ_func_0_0(0) * occ_func_0_0(48) + occ_func_0_0(0) * occ_func_0_0(44) * occ_func_0_0(43) + occ_func_0_0(0) * occ_func_0_0(53) * occ_func_0_0(50) + occ_func_0_0(47) * occ_func_0_0(52) * occ_func_0_0(0) + occ_func_0_0(44) * occ_func_0_0(0) * occ_func_0_0(45) + occ_func_0_0(54) * occ_func_0_0(0) * occ_func_0_0(51) + occ_func_0_0(52) * occ_func_0_0(46) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(43) * occ_func_0_0(45) + occ_func_0_0(51) * occ_func_0_0(0) * occ_func_0_0(50) + occ_func_0_0(49) * occ_func_0_0(47) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(46) * occ_func_0_0(48)) / 8.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_34_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(49) * occ_func_0_0(51) + occ_func_0_0(48) * occ_func_0_0(50) + occ_func_0_0(46) * occ_func_0_0(47) + occ_func_0_0(52) * occ_func_0_0(54) + occ_func_0_0(51) * occ_func_0_0(45) + occ_func_0_0(46) * occ_func_0_0(43) + occ_func_0_0(53) * occ_func_0_0(52) + occ_func_0_0(50) * occ_func_0_0(45) + occ_func_0_0(47) * occ_func_0_0(44) + occ_func_0_0(53) * occ_func_0_0(54) + occ_func_0_0(44) * occ_func_0_0(49) + occ_func_0_0(48) * occ_func_0_0(43) + occ_func_0_0(54) * occ_func_0_0(49) + occ_func_0_0(53) * occ_func_0_0(48) + occ_func_0_0(44) * occ_func_0_0(43) + occ_func_0_0(53) * occ_func_0_0(50) + occ_func_0_0(47) * occ_func_0_0(52) + occ_func_0_0(44) * occ_func_0_0(45) + occ_func_0_0(54) * occ_func_0_0(51) + occ_func_0_0(52) * occ_func_0_0(46) + occ_func_0_0(43) * occ_func_0_0(45) + occ_func_0_0(51) * occ_func_0_0(50) + occ_func_0_0(49) * occ_func_0_0(47) + occ_func_0_0(46) * occ_func_0_0(48)) / 8.;
  }

  /**** Basis functions for orbit 35****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 0.0000000 1.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_35_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_0(12) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(11)) / 2.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_35_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6)) / 2.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_35_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(11) + occ_func_0_0(9) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(6)) / 2.;
  }

  /**** Basis functions for orbit 36****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 -1.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_36_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(9) * occ_func_0_0(12) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(18) + occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(12) * occ_func_0_0(11) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(9) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(7) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(12) * occ_func_0_0(0)) / 12.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_36_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(9) * occ_func_0_0(12) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(18) + occ_func_0_0(10) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(12) * occ_func_0_0(11) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(16) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(7) * occ_func_0_0(5) + occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(13) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_0(10) + occ_func_0_0(3) * occ_func_0_0(14) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(9) + occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(4) * occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(15) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(7) * occ_func_0_0(15) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_0(6) + occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_0(9) + occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(13) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(3) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(13) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(15) + occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_0(3) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(14) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(1)) / 12.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_36_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(8) * occ_func_0_0(9) * occ_func_0_0(12) + occ_func_0_0(5) * occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(12) * occ_func_0_0(9) * occ_func_0_0(18) + occ_func_0_0(10) * occ_func_0_0(4) * occ_func_0_0(11) + occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(2) + occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(12) + occ_func_0_0(10) * occ_func_0_0(5) * occ_func_0_0(17) + occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(1) + occ_func_0_0(9) * occ_func_0_0(12) * occ_func_0_0(11) + occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(16) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(7) + occ_func_0_0(1) * occ_func_0_0(3) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(9) + occ_func_0_0(12) * occ_func_0_0(7) * occ_func_0_0(5) + occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(13) + occ_func_0_0(7) * occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_0(10) + occ_func_0_0(3) * occ_func_0_0(14) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(1) * occ_func_0_0(5) + occ_func_0_0(11) * occ_func_0_0(8) * occ_func_0_0(9) + occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_0(7) + occ_func_0_0(4) * occ_func_0_0(16) * occ_func_0_0(6) + occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(15) + occ_func_0_0(10) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(7) * occ_func_0_0(15) * occ_func_0_0(9) + occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_0(6) + occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(4) + occ_func_0_0(11) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_0(9) + occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(13) * occ_func_0_0(5) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(8) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(5) + occ_func_0_0(9) * occ_func_0_0(3) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(13) + occ_func_0_0(10) * occ_func_0_0(12) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(11) + occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(15) + occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_0(2) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(7) * occ_func_0_0(9) * occ_func_0_0(3) + occ_func_0_0(6) * occ_func_0_0(8) * occ_func_0_0(14) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(1)) / 12.;
  }

  /**** Basis functions for orbit 37****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 -1.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

1.0000000 1.0000000 -1.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::eval_bfunc_37_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(11) * occ_func_0_0(18)) / 3.;
  }

  template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_eval_bfunc_37_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(13) * occ_func_0_0(1) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) + occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(15) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(11) * occ_func_0_0(18) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(13) * occ_func_0_0(2) * occ_func_0_0(4) * occ_func_0_0(0)) / 3.;
  }

    template<typename Scalar>
  Scalar fcc_binary_Clexulator_default::site_deval_bfunc_37_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(8) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(5) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(8) + occ_func_0_0(13) * occ_func_0_0(1) * occ_func_0_0(5) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(7) + occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(10) + occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(15) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(9) * occ_func_0_0(11) * occ_func_0_0(18) + occ_func_0_0(4) * occ_func_0_0(16) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(9) + occ_func_0_0(13) * occ_func_0_0(2) * occ_func_0_0(4)) / 3.;
  }

} // namespace clexulator
} // namespace CASM


extern "C" {
  /// \brief Returns a clexulator::BaseClexulator* owning a fcc_binary_Clexulator_default
  CASM::clexulator::BaseClexulator *make_fcc_binary_Clexulator_default() {
    return new CASM::clexulator::fcc_binary_Clexulator_default();
  }

}

