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
    [ 1.000000000000, 0.000000000000, 0.000000000000 ],
    [ 0.000000000000, 1.000000000000, 0.000000000000 ],
    [ 0.000000000000, 0.000000000000, 10.000000000000 ]
  ],
  "title" : "square_binary"
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
      "generating_group" : [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 ],
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
          "max_length" : 3.010000000000
        },
        "4" : {
          "max_length" : 2.510000000000
        }
      }
    }
  }
}

**/


/// \brief Returns a clexulator::BaseClexulator* owning a square_binary_Clexulator_default
extern "C" CASM::clexulator::BaseClexulator *make_square_binary_Clexulator_default();

namespace CASM {
namespace clexulator {

  /****** GENERATED CLEXPARAMPACK DEFINITION ******/


  typedef BasicClexParamPack ParamPack;


  /****** GENERATED CLEXULATOR DEFINITION ******/

  class square_binary_Clexulator_default : public clexulator::BaseClexulator {

  public:

    square_binary_Clexulator_default();

    ~square_binary_Clexulator_default();

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


  private:

    // ParamPack object, which stores temporary data for calculations
    mutable ParamPack m_params;

    // typedef for method pointers of scalar type double
    typedef double (square_binary_Clexulator_default::*BasisFuncPtr_0)() const;

    // typedef for method pointers
    typedef double (square_binary_Clexulator_default::*DeltaBasisFuncPtr_0)(int, int) const;

    // array of pointers to member functions for calculating basis functions of scalar type double
    BasisFuncPtr_0 m_orbit_func_table_0[33];

    // array of pointers to member functions for calculating flower functions of scalar type double
    BasisFuncPtr_0 m_flower_func_table_0[1][33];

    // array of pointers to member functions for calculating DELTA flower functions of scalar type double
    DeltaBasisFuncPtr_0 m_delta_func_table_0[1][33];

    // Occupation Function tables for basis sites in asymmetric unit 0:
    //   - basis site 0:
    double m_occ_func_0_0[2];

    //ClexParamPack allocation for evaluated correlations 
    ParamPack::Key m_corr_param_key;
    //ClexParamPack allocation for DoF occ
    ParamPack::Key m_occ_site_func_param_key;

    /// \brief Clone the square_binary_Clexulator_default
    BaseClexulator *_clone() const override {
      return new square_binary_Clexulator_default(*this);
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

  square_binary_Clexulator_default::square_binary_Clexulator_default() :
    BaseClexulator(29, 33, 1) {
    m_occ_func_0_0[0] = -1.0000000000, m_occ_func_0_0[1] = 1.0000000000;

    m_occ_site_func_param_key = m_params.allocate("occ_site_func", 1, 29, true);

    m_corr_param_key = m_params.allocate("corr", corr_size(), 1, false);

    m_orbit_func_table_0[0] = &square_binary_Clexulator_default::eval_bfunc_0_0<double>;
    m_orbit_func_table_0[1] = &square_binary_Clexulator_default::eval_bfunc_1_0<double>;
    m_orbit_func_table_0[2] = &square_binary_Clexulator_default::eval_bfunc_2_0<double>;
    m_orbit_func_table_0[3] = &square_binary_Clexulator_default::eval_bfunc_3_0<double>;
    m_orbit_func_table_0[4] = &square_binary_Clexulator_default::eval_bfunc_4_0<double>;
    m_orbit_func_table_0[5] = &square_binary_Clexulator_default::eval_bfunc_5_0<double>;
    m_orbit_func_table_0[6] = &square_binary_Clexulator_default::eval_bfunc_6_0<double>;
    m_orbit_func_table_0[7] = &square_binary_Clexulator_default::eval_bfunc_7_0<double>;
    m_orbit_func_table_0[8] = &square_binary_Clexulator_default::eval_bfunc_8_0<double>;
    m_orbit_func_table_0[9] = &square_binary_Clexulator_default::eval_bfunc_9_0<double>;
    m_orbit_func_table_0[10] = &square_binary_Clexulator_default::eval_bfunc_10_0<double>;
    m_orbit_func_table_0[11] = &square_binary_Clexulator_default::eval_bfunc_11_0<double>;
    m_orbit_func_table_0[12] = &square_binary_Clexulator_default::eval_bfunc_12_0<double>;
    m_orbit_func_table_0[13] = &square_binary_Clexulator_default::eval_bfunc_13_0<double>;
    m_orbit_func_table_0[14] = &square_binary_Clexulator_default::eval_bfunc_14_0<double>;
    m_orbit_func_table_0[15] = &square_binary_Clexulator_default::eval_bfunc_15_0<double>;
    m_orbit_func_table_0[16] = &square_binary_Clexulator_default::eval_bfunc_16_0<double>;
    m_orbit_func_table_0[17] = &square_binary_Clexulator_default::eval_bfunc_17_0<double>;
    m_orbit_func_table_0[18] = &square_binary_Clexulator_default::eval_bfunc_18_0<double>;
    m_orbit_func_table_0[19] = &square_binary_Clexulator_default::eval_bfunc_19_0<double>;
    m_orbit_func_table_0[20] = &square_binary_Clexulator_default::eval_bfunc_20_0<double>;
    m_orbit_func_table_0[21] = &square_binary_Clexulator_default::eval_bfunc_21_0<double>;
    m_orbit_func_table_0[22] = &square_binary_Clexulator_default::eval_bfunc_22_0<double>;
    m_orbit_func_table_0[23] = &square_binary_Clexulator_default::eval_bfunc_23_0<double>;
    m_orbit_func_table_0[24] = &square_binary_Clexulator_default::eval_bfunc_24_0<double>;
    m_orbit_func_table_0[25] = &square_binary_Clexulator_default::eval_bfunc_25_0<double>;
    m_orbit_func_table_0[26] = &square_binary_Clexulator_default::eval_bfunc_26_0<double>;
    m_orbit_func_table_0[27] = &square_binary_Clexulator_default::eval_bfunc_27_0<double>;
    m_orbit_func_table_0[28] = &square_binary_Clexulator_default::eval_bfunc_28_0<double>;
    m_orbit_func_table_0[29] = &square_binary_Clexulator_default::eval_bfunc_29_0<double>;
    m_orbit_func_table_0[30] = &square_binary_Clexulator_default::eval_bfunc_30_0<double>;
    m_orbit_func_table_0[31] = &square_binary_Clexulator_default::eval_bfunc_31_0<double>;
    m_orbit_func_table_0[32] = &square_binary_Clexulator_default::eval_bfunc_32_0<double>;


    m_flower_func_table_0[0][0] = &square_binary_Clexulator_default::zero_func<double>;
    m_flower_func_table_0[0][1] = &square_binary_Clexulator_default::site_eval_bfunc_1_0_at_0<double>;
    m_flower_func_table_0[0][2] = &square_binary_Clexulator_default::site_eval_bfunc_2_0_at_0<double>;
    m_flower_func_table_0[0][3] = &square_binary_Clexulator_default::site_eval_bfunc_3_0_at_0<double>;
    m_flower_func_table_0[0][4] = &square_binary_Clexulator_default::site_eval_bfunc_4_0_at_0<double>;
    m_flower_func_table_0[0][5] = &square_binary_Clexulator_default::site_eval_bfunc_5_0_at_0<double>;
    m_flower_func_table_0[0][6] = &square_binary_Clexulator_default::site_eval_bfunc_6_0_at_0<double>;
    m_flower_func_table_0[0][7] = &square_binary_Clexulator_default::site_eval_bfunc_7_0_at_0<double>;
    m_flower_func_table_0[0][8] = &square_binary_Clexulator_default::site_eval_bfunc_8_0_at_0<double>;
    m_flower_func_table_0[0][9] = &square_binary_Clexulator_default::site_eval_bfunc_9_0_at_0<double>;
    m_flower_func_table_0[0][10] = &square_binary_Clexulator_default::site_eval_bfunc_10_0_at_0<double>;
    m_flower_func_table_0[0][11] = &square_binary_Clexulator_default::site_eval_bfunc_11_0_at_0<double>;
    m_flower_func_table_0[0][12] = &square_binary_Clexulator_default::site_eval_bfunc_12_0_at_0<double>;
    m_flower_func_table_0[0][13] = &square_binary_Clexulator_default::site_eval_bfunc_13_0_at_0<double>;
    m_flower_func_table_0[0][14] = &square_binary_Clexulator_default::site_eval_bfunc_14_0_at_0<double>;
    m_flower_func_table_0[0][15] = &square_binary_Clexulator_default::site_eval_bfunc_15_0_at_0<double>;
    m_flower_func_table_0[0][16] = &square_binary_Clexulator_default::site_eval_bfunc_16_0_at_0<double>;
    m_flower_func_table_0[0][17] = &square_binary_Clexulator_default::site_eval_bfunc_17_0_at_0<double>;
    m_flower_func_table_0[0][18] = &square_binary_Clexulator_default::site_eval_bfunc_18_0_at_0<double>;
    m_flower_func_table_0[0][19] = &square_binary_Clexulator_default::site_eval_bfunc_19_0_at_0<double>;
    m_flower_func_table_0[0][20] = &square_binary_Clexulator_default::site_eval_bfunc_20_0_at_0<double>;
    m_flower_func_table_0[0][21] = &square_binary_Clexulator_default::site_eval_bfunc_21_0_at_0<double>;
    m_flower_func_table_0[0][22] = &square_binary_Clexulator_default::site_eval_bfunc_22_0_at_0<double>;
    m_flower_func_table_0[0][23] = &square_binary_Clexulator_default::site_eval_bfunc_23_0_at_0<double>;
    m_flower_func_table_0[0][24] = &square_binary_Clexulator_default::site_eval_bfunc_24_0_at_0<double>;
    m_flower_func_table_0[0][25] = &square_binary_Clexulator_default::site_eval_bfunc_25_0_at_0<double>;
    m_flower_func_table_0[0][26] = &square_binary_Clexulator_default::site_eval_bfunc_26_0_at_0<double>;
    m_flower_func_table_0[0][27] = &square_binary_Clexulator_default::site_eval_bfunc_27_0_at_0<double>;
    m_flower_func_table_0[0][28] = &square_binary_Clexulator_default::site_eval_bfunc_28_0_at_0<double>;
    m_flower_func_table_0[0][29] = &square_binary_Clexulator_default::site_eval_bfunc_29_0_at_0<double>;
    m_flower_func_table_0[0][30] = &square_binary_Clexulator_default::site_eval_bfunc_30_0_at_0<double>;
    m_flower_func_table_0[0][31] = &square_binary_Clexulator_default::site_eval_bfunc_31_0_at_0<double>;
    m_flower_func_table_0[0][32] = &square_binary_Clexulator_default::site_eval_bfunc_32_0_at_0<double>;


    m_delta_func_table_0[0][0] = &square_binary_Clexulator_default::zero_func<double>;
    m_delta_func_table_0[0][1] = &square_binary_Clexulator_default::site_deval_bfunc_1_0_at_0<double>;
    m_delta_func_table_0[0][2] = &square_binary_Clexulator_default::site_deval_bfunc_2_0_at_0<double>;
    m_delta_func_table_0[0][3] = &square_binary_Clexulator_default::site_deval_bfunc_3_0_at_0<double>;
    m_delta_func_table_0[0][4] = &square_binary_Clexulator_default::site_deval_bfunc_4_0_at_0<double>;
    m_delta_func_table_0[0][5] = &square_binary_Clexulator_default::site_deval_bfunc_5_0_at_0<double>;
    m_delta_func_table_0[0][6] = &square_binary_Clexulator_default::site_deval_bfunc_6_0_at_0<double>;
    m_delta_func_table_0[0][7] = &square_binary_Clexulator_default::site_deval_bfunc_7_0_at_0<double>;
    m_delta_func_table_0[0][8] = &square_binary_Clexulator_default::site_deval_bfunc_8_0_at_0<double>;
    m_delta_func_table_0[0][9] = &square_binary_Clexulator_default::site_deval_bfunc_9_0_at_0<double>;
    m_delta_func_table_0[0][10] = &square_binary_Clexulator_default::site_deval_bfunc_10_0_at_0<double>;
    m_delta_func_table_0[0][11] = &square_binary_Clexulator_default::site_deval_bfunc_11_0_at_0<double>;
    m_delta_func_table_0[0][12] = &square_binary_Clexulator_default::site_deval_bfunc_12_0_at_0<double>;
    m_delta_func_table_0[0][13] = &square_binary_Clexulator_default::site_deval_bfunc_13_0_at_0<double>;
    m_delta_func_table_0[0][14] = &square_binary_Clexulator_default::site_deval_bfunc_14_0_at_0<double>;
    m_delta_func_table_0[0][15] = &square_binary_Clexulator_default::site_deval_bfunc_15_0_at_0<double>;
    m_delta_func_table_0[0][16] = &square_binary_Clexulator_default::site_deval_bfunc_16_0_at_0<double>;
    m_delta_func_table_0[0][17] = &square_binary_Clexulator_default::site_deval_bfunc_17_0_at_0<double>;
    m_delta_func_table_0[0][18] = &square_binary_Clexulator_default::site_deval_bfunc_18_0_at_0<double>;
    m_delta_func_table_0[0][19] = &square_binary_Clexulator_default::site_deval_bfunc_19_0_at_0<double>;
    m_delta_func_table_0[0][20] = &square_binary_Clexulator_default::site_deval_bfunc_20_0_at_0<double>;
    m_delta_func_table_0[0][21] = &square_binary_Clexulator_default::site_deval_bfunc_21_0_at_0<double>;
    m_delta_func_table_0[0][22] = &square_binary_Clexulator_default::site_deval_bfunc_22_0_at_0<double>;
    m_delta_func_table_0[0][23] = &square_binary_Clexulator_default::site_deval_bfunc_23_0_at_0<double>;
    m_delta_func_table_0[0][24] = &square_binary_Clexulator_default::site_deval_bfunc_24_0_at_0<double>;
    m_delta_func_table_0[0][25] = &square_binary_Clexulator_default::site_deval_bfunc_25_0_at_0<double>;
    m_delta_func_table_0[0][26] = &square_binary_Clexulator_default::site_deval_bfunc_26_0_at_0<double>;
    m_delta_func_table_0[0][27] = &square_binary_Clexulator_default::site_deval_bfunc_27_0_at_0<double>;
    m_delta_func_table_0[0][28] = &square_binary_Clexulator_default::site_deval_bfunc_28_0_at_0<double>;
    m_delta_func_table_0[0][29] = &square_binary_Clexulator_default::site_deval_bfunc_29_0_at_0<double>;
    m_delta_func_table_0[0][30] = &square_binary_Clexulator_default::site_deval_bfunc_30_0_at_0<double>;
    m_delta_func_table_0[0][31] = &square_binary_Clexulator_default::site_deval_bfunc_31_0_at_0<double>;
    m_delta_func_table_0[0][32] = &square_binary_Clexulator_default::site_deval_bfunc_32_0_at_0<double>;


    m_weight_matrix.row(0) << 1, 0, 0;
    m_weight_matrix.row(1) << 0, 1, 0;
    m_weight_matrix.row(2) << 0, 0, 100;

    m_sublat_indices = std::set<int>{0};

    m_n_sublattices = 1;

    m_neighborhood = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0),
      xtal::UnitCell(2, 2, 0),
      xtal::UnitCell(3, 0, 0)
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
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, 0, 0)
    };

    m_orbit_site_neighborhood[2] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0)
    };

    m_orbit_neighborhood[3] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 1, 0)
    };

    m_orbit_site_neighborhood[3] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 1, 0)
    };

    m_orbit_neighborhood[4] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[4] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[5] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[5] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[6] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 2, 0)
    };

    m_orbit_site_neighborhood[6] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -2, 0),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 2, 0)
    };

    m_orbit_neighborhood[7] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(3, 0, 0)
    };

    m_orbit_site_neighborhood[7] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 0, 0),
      xtal::UnitCellCoord(0, 0, -3, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 3, 0),
      xtal::UnitCellCoord(0, 3, 0, 0)
    };

    m_orbit_neighborhood[8] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0)
    };

    m_orbit_site_neighborhood[8] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0)
    };

    m_orbit_neighborhood[9] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[9] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[10] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[10] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[11] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[11] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[12] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[12] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[13] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[13] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[14] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[14] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[15] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 2, 0)
    };

    m_orbit_site_neighborhood[15] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -2, 0),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 2, 0)
    };

    m_orbit_neighborhood[16] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 2, 0)
    };

    m_orbit_site_neighborhood[16] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -2, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 2, 0)
    };

    m_orbit_neighborhood[17] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0),
      xtal::UnitCell(2, 2, 0)
    };

    m_orbit_site_neighborhood[17] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -2, 0),
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0),
      xtal::UnitCellCoord(0, 2, 2, 0)
    };

    m_orbit_neighborhood[18] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(3, 0, 0)
    };

    m_orbit_site_neighborhood[18] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 0, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, 0, -3, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 0, 3, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 3, 0, 0)
    };

    m_orbit_neighborhood[19] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0),
      xtal::UnitCell(3, 0, 0)
    };

    m_orbit_site_neighborhood[19] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 0, 0),
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -3, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 3, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0),
      xtal::UnitCellCoord(0, 3, 0, 0)
    };

    m_orbit_neighborhood[20] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-3, 0, 0),
      xtal::UnitCell(-2, -2, 0),
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-2, 2, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -3, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 3, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0),
      xtal::UnitCell(2, 2, 0),
      xtal::UnitCell(3, 0, 0)
    };

    m_orbit_site_neighborhood[20] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -3, 0, 0),
      xtal::UnitCellCoord(0, -2, -2, 0),
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -2, 2, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -3, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 3, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0),
      xtal::UnitCellCoord(0, 2, 2, 0),
      xtal::UnitCellCoord(0, 3, 0, 0)
    };

    m_orbit_neighborhood[21] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0)
    };

    m_orbit_site_neighborhood[21] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0)
    };

    m_orbit_neighborhood[22] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[22] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[23] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(2, 0, 0)
    };

    m_orbit_site_neighborhood[23] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0)
    };

    m_orbit_neighborhood[24] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[24] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[25] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[25] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[26] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[26] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[27] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[27] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[28] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[28] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[29] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[29] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[30] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[30] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[31] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[31] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

    m_orbit_neighborhood[32] = std::set<xtal::UnitCell> {
      xtal::UnitCell(-2, -1, 0),
      xtal::UnitCell(-2, 0, 0),
      xtal::UnitCell(-2, 1, 0),
      xtal::UnitCell(-1, -2, 0),
      xtal::UnitCell(-1, -1, 0),
      xtal::UnitCell(-1, 0, 0),
      xtal::UnitCell(-1, 1, 0),
      xtal::UnitCell(-1, 2, 0),
      xtal::UnitCell(0, -2, 0),
      xtal::UnitCell(0, -1, 0),
      xtal::UnitCell(0, 0, 0),
      xtal::UnitCell(0, 1, 0),
      xtal::UnitCell(0, 2, 0),
      xtal::UnitCell(1, -2, 0),
      xtal::UnitCell(1, -1, 0),
      xtal::UnitCell(1, 0, 0),
      xtal::UnitCell(1, 1, 0),
      xtal::UnitCell(1, 2, 0),
      xtal::UnitCell(2, -1, 0),
      xtal::UnitCell(2, 0, 0),
      xtal::UnitCell(2, 1, 0)
    };

    m_orbit_site_neighborhood[32] = std::set<xtal::UnitCellCoord> {
      xtal::UnitCellCoord(0, -2, -1, 0),
      xtal::UnitCellCoord(0, -2, 0, 0),
      xtal::UnitCellCoord(0, -2, 1, 0),
      xtal::UnitCellCoord(0, -1, -2, 0),
      xtal::UnitCellCoord(0, -1, -1, 0),
      xtal::UnitCellCoord(0, -1, 0, 0),
      xtal::UnitCellCoord(0, -1, 1, 0),
      xtal::UnitCellCoord(0, -1, 2, 0),
      xtal::UnitCellCoord(0, 0, -2, 0),
      xtal::UnitCellCoord(0, 0, -1, 0),
      xtal::UnitCellCoord(0, 0, 0, 0),
      xtal::UnitCellCoord(0, 0, 1, 0),
      xtal::UnitCellCoord(0, 0, 2, 0),
      xtal::UnitCellCoord(0, 1, -2, 0),
      xtal::UnitCellCoord(0, 1, -1, 0),
      xtal::UnitCellCoord(0, 1, 0, 0),
      xtal::UnitCellCoord(0, 1, 1, 0),
      xtal::UnitCellCoord(0, 1, 2, 0),
      xtal::UnitCellCoord(0, 2, -1, 0),
      xtal::UnitCellCoord(0, 2, 0, 0),
      xtal::UnitCellCoord(0, 2, 1, 0)
    };

  }


  square_binary_Clexulator_default::~square_binary_Clexulator_default() {
    //nothing here for now
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void square_binary_Clexulator_default::_calc_global_corr_contribution(double *corr_begin) const {
    _calc_global_corr_contribution();
    for(size_type i = 0; i < corr_size(); i++) {
      *(corr_begin + i) = ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
    }
  }

  /// \brief Calculate contribution to global correlations from one unit cell
  void square_binary_Clexulator_default::_calc_global_corr_contribution() const {
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
  void square_binary_Clexulator_default::_calc_restricted_global_corr_contribution(double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    _calc_restricted_global_corr_contribution(ind_list_begin, ind_list_end);
    for(; ind_list_begin < ind_list_end; ind_list_begin++) {
      *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(m_params, m_corr_param_key, *ind_list_begin);
    }
  }

  /// \brief Calculate contribution to select global correlations from one unit cell
  void square_binary_Clexulator_default::_calc_restricted_global_corr_contribution(size_type const *ind_list_begin, size_type const *ind_list_end) const {
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
  void square_binary_Clexulator_default::_calc_point_corr(int nlist_ind, double *corr_begin) const {
    _calc_point_corr(nlist_ind);
    for(size_type i = 0; i < corr_size(); i++) {
      *(corr_begin + i) = ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
    }
  }

  /// \brief Calculate point correlations about basis site 'nlist_ind'
  void square_binary_Clexulator_default::_calc_point_corr(int nlist_ind) const {
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
  void square_binary_Clexulator_default::_calc_restricted_point_corr(int nlist_ind, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    _calc_restricted_point_corr(nlist_ind, ind_list_begin, ind_list_end);
    for(; ind_list_begin < ind_list_end; ind_list_begin++) {
      *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(m_params, m_corr_param_key, *ind_list_begin);
    }
  }

  /// \brief Calculate select point correlations about basis site 'nlist_ind'
  void square_binary_Clexulator_default::_calc_restricted_point_corr(int nlist_ind, size_type const *ind_list_begin, size_type const *ind_list_end) const {
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
  void square_binary_Clexulator_default::_calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double *corr_begin) const {
    _calc_delta_point_corr(nlist_ind, occ_i, occ_f);
    for(size_type i = 0; i < corr_size(); i++) {
      *(corr_begin + i) = ParamPack::Val<double>::get(m_params, m_corr_param_key, i);
    }
  }

  /// \brief Calculate the change in point correlations due to changing an occupant
  void square_binary_Clexulator_default::_calc_delta_point_corr(int nlist_ind, int occ_i, int occ_f) const {
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
  void square_binary_Clexulator_default::_calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f, double *corr_begin, size_type const *ind_list_begin, size_type const *ind_list_end) const {
    _calc_restricted_delta_point_corr(nlist_ind, occ_i, occ_f, ind_list_begin, ind_list_end);
    for(; ind_list_begin < ind_list_end; ind_list_begin++) {
      *(corr_begin + *ind_list_begin) = ParamPack::Val<double>::get(m_params, m_corr_param_key, *ind_list_begin);
    }
  }

  /// \brief Calculate the change in select point correlations due to changing an occupant
  void square_binary_Clexulator_default::_calc_restricted_delta_point_corr(int nlist_ind, int occ_i, int occ_f, size_type const *ind_list_begin, size_type const *ind_list_end) const {
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
  void square_binary_Clexulator_default::_point_prepare(int nlist_ind) const {
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
      }
      break;
    }
  }
  template<typename Scalar>
  void square_binary_Clexulator_default::_global_prepare() const {
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
  }
  }

  // Basis functions for empty cluster:
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_0_0() const {
    return 1;
  }

  /**** Basis functions for orbit 1****
0.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_1_0() const {
    return occ_func_0_0(0);
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_1_0_at_0() const {
    return occ_func_0_0(0);
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_1_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]);
  }

  /**** Basis functions for orbit 2****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_2_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(4) * occ_func_0_0(0)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_2_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_2_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) + occ_func_0_0(2) + occ_func_0_0(4) + occ_func_0_0(1)) / 2.;
  }

  /**** Basis functions for orbit 3****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_3_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(0) * occ_func_0_0(8)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_3_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(0)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_3_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) + occ_func_0_0(6) + occ_func_0_0(8) + occ_func_0_0(5)) / 2.;
  }

  /**** Basis functions for orbit 4****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_4_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(12) * occ_func_0_0(0)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_4_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_4_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(11) + occ_func_0_0(10) + occ_func_0_0(12) + occ_func_0_0(9)) / 2.;
  }

  /**** Basis functions for orbit 5****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_5_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(19)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_5_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(14) * occ_func_0_0(0)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_5_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(17) + occ_func_0_0(16) + occ_func_0_0(20) + occ_func_0_0(13) + occ_func_0_0(18) + occ_func_0_0(15) + occ_func_0_0(19) + occ_func_0_0(14)) / 4.;
  }

  /**** Basis functions for orbit 6****
0.0000000 0.0000000 0.0000000 A  B  

2.0000000 -2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_6_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(24)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_6_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(21) * occ_func_0_0(0)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_6_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(23) + occ_func_0_0(22) + occ_func_0_0(24) + occ_func_0_0(21)) / 2.;
  }

  /**** Basis functions for orbit 7****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 3.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_7_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(27) + occ_func_0_0(28) * occ_func_0_0(0)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_7_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(27) + occ_func_0_0(26) * occ_func_0_0(0) + occ_func_0_0(28) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_7_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(27) + occ_func_0_0(26) + occ_func_0_0(28) + occ_func_0_0(25)) / 2.;
  }

  /**** Basis functions for orbit 8****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_8_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_8_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_8_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(4) + occ_func_0_0(5) * occ_func_0_0(1) + occ_func_0_0(4) * occ_func_0_0(3) + occ_func_0_0(7) * occ_func_0_0(2) + occ_func_0_0(1) * occ_func_0_0(6) + occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(3) + occ_func_0_0(1) * occ_func_0_0(2) + occ_func_0_0(8) * occ_func_0_0(4) + occ_func_0_0(3) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(5)) / 4.;
  }

  /**** Basis functions for orbit 9****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_9_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(0)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_9_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(9)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_9_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(2) + occ_func_0_0(12) * occ_func_0_0(4) + occ_func_0_0(4) * occ_func_0_0(1) + occ_func_0_0(1) * occ_func_0_0(9)) / 2.;
  }

  /**** Basis functions for orbit 10****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_10_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(12) * occ_func_0_0(8) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_10_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(9) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(9)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_10_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(8) + occ_func_0_0(9) * occ_func_0_0(5) + occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(7) + occ_func_0_0(5) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(6) + occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(8) + occ_func_0_0(7) * occ_func_0_0(5) + occ_func_0_0(6) * occ_func_0_0(9)) / 4.;
  }

  /**** Basis functions for orbit 11****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_11_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_11_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(5) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(14)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_11_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(7) + occ_func_0_0(2) * occ_func_0_0(17) + occ_func_0_0(6) * occ_func_0_0(16) + occ_func_0_0(4) * occ_func_0_0(20) + occ_func_0_0(1) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(13) + occ_func_0_0(8) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(5) + occ_func_0_0(1) * occ_func_0_0(13) + occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(3) * occ_func_0_0(16) + occ_func_0_0(2) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(3) * occ_func_0_0(5) + occ_func_0_0(2) * occ_func_0_0(15) + occ_func_0_0(3) * occ_func_0_0(18) + occ_func_0_0(2) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(15) + occ_func_0_0(4) * occ_func_0_0(19) + occ_func_0_0(1) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(14) + occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(4) * occ_func_0_0(6) + occ_func_0_0(1) * occ_func_0_0(14)) / 8.;
  }

  /**** Basis functions for orbit 12****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_12_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(12) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(11)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_12_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(12) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(9) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(12) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(14) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(9) + occ_func_0_0(12) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(9) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(12) + occ_func_0_0(9) * occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(11) * occ_func_0_0(16) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(10) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(16)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_12_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(12) + occ_func_0_0(2) * occ_func_0_0(19) + occ_func_0_0(9) * occ_func_0_0(14) + occ_func_0_0(4) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(15) + occ_func_0_0(11) * occ_func_0_0(18) + occ_func_0_0(4) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(15) + occ_func_0_0(12) * occ_func_0_0(19) + occ_func_0_0(3) * occ_func_0_0(14) + occ_func_0_0(2) * occ_func_0_0(9) + occ_func_0_0(12) * occ_func_0_0(20) + occ_func_0_0(3) * occ_func_0_0(9) + occ_func_0_0(2) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(20) + occ_func_0_0(2) * occ_func_0_0(12) + occ_func_0_0(9) * occ_func_0_0(13) + occ_func_0_0(4) * occ_func_0_0(17) + occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(1) * occ_func_0_0(10) + occ_func_0_0(4) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(17) + occ_func_0_0(1) * occ_func_0_0(16)) / 8.;
  }

  /**** Basis functions for orbit 13****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 0.0000000 A  B  

2.0000000 -1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_13_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(18) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_13_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(19) + occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(14) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(18) + occ_func_0_0(15) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(18) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(15) + occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(14)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_13_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(17) * occ_func_0_0(19) + occ_func_0_0(16) * occ_func_0_0(8) + occ_func_0_0(14) * occ_func_0_0(5) + occ_func_0_0(20) * occ_func_0_0(18) + occ_func_0_0(15) * occ_func_0_0(7) + occ_func_0_0(13) * occ_func_0_0(6) + occ_func_0_0(20) * occ_func_0_0(7) + occ_func_0_0(18) * occ_func_0_0(6) + occ_func_0_0(13) * occ_func_0_0(15) + occ_func_0_0(19) * occ_func_0_0(8) + occ_func_0_0(17) * occ_func_0_0(5) + occ_func_0_0(16) * occ_func_0_0(14)) / 4.;
  }

  /**** Basis functions for orbit 14****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_14_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(11) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_14_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(9) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(11) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(9)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_14_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(16) * occ_func_0_0(18) + occ_func_0_0(9) * occ_func_0_0(15) + occ_func_0_0(20) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(19) + occ_func_0_0(13) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(19) + occ_func_0_0(11) * occ_func_0_0(14) + occ_func_0_0(13) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(17) * occ_func_0_0(15) + occ_func_0_0(16) * occ_func_0_0(9)) / 4.;
  }

  /**** Basis functions for orbit 15****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

2.0000000 -2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_15_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(24)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_15_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(23) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(22) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(24) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(21) * occ_func_0_0(5) * occ_func_0_0(0)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_15_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(23) + occ_func_0_0(6) * occ_func_0_0(7) + occ_func_0_0(22) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(24) + occ_func_0_0(5) * occ_func_0_0(8) + occ_func_0_0(21) * occ_func_0_0(5)) / 2.;
  }

  /**** Basis functions for orbit 16****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

2.0000000 2.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_16_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(24) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(23) + occ_func_0_0(24) * occ_func_0_0(12) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_16_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(24) + occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(21) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(23) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(22) + occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(23) + occ_func_0_0(22) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(24) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(21)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_16_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(11) * occ_func_0_0(24) + occ_func_0_0(10) * occ_func_0_0(12) + occ_func_0_0(21) * occ_func_0_0(9) + occ_func_0_0(12) * occ_func_0_0(11) + occ_func_0_0(23) * occ_func_0_0(10) + occ_func_0_0(9) * occ_func_0_0(22) + occ_func_0_0(12) * occ_func_0_0(23) + occ_func_0_0(22) * occ_func_0_0(11) + occ_func_0_0(9) * occ_func_0_0(10) + occ_func_0_0(24) * occ_func_0_0(12) + occ_func_0_0(11) * occ_func_0_0(9) + occ_func_0_0(10) * occ_func_0_0(21)) / 4.;
  }

  /**** Basis functions for orbit 17****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

2.0000000 -1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_17_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(19) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(18) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(20) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(17) * occ_func_0_0(23) * occ_func_0_0(0)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_17_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(19) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(14) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(18) + occ_func_0_0(15) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(15) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(19) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(22) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(13) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(20) + occ_func_0_0(13) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(17) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(16) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(22)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_17_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(19) + occ_func_0_0(2) * occ_func_0_0(23) + occ_func_0_0(14) * occ_func_0_0(22) + occ_func_0_0(4) * occ_func_0_0(24) + occ_func_0_0(1) * occ_func_0_0(18) + occ_func_0_0(15) * occ_func_0_0(21) + occ_func_0_0(18) * occ_func_0_0(24) + occ_func_0_0(4) * occ_func_0_0(15) + occ_func_0_0(1) * occ_func_0_0(21) + occ_func_0_0(19) * occ_func_0_0(23) + occ_func_0_0(3) * occ_func_0_0(22) + occ_func_0_0(2) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(24) + occ_func_0_0(3) * occ_func_0_0(13) + occ_func_0_0(2) * occ_func_0_0(21) + occ_func_0_0(3) * occ_func_0_0(24) + occ_func_0_0(2) * occ_func_0_0(20) + occ_func_0_0(13) * occ_func_0_0(21) + occ_func_0_0(4) * occ_func_0_0(23) + occ_func_0_0(1) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(22) + occ_func_0_0(17) * occ_func_0_0(23) + occ_func_0_0(4) * occ_func_0_0(16) + occ_func_0_0(1) * occ_func_0_0(22)) / 8.;
  }

  /**** Basis functions for orbit 18****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

0.0000000 3.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_18_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(27) + occ_func_0_0(28) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(28) + occ_func_0_0(27) * occ_func_0_0(11) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_18_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(27) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(26) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(28) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(25) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(28) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(25) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(27) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(26)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_18_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(27) + occ_func_0_0(2) * occ_func_0_0(11) + occ_func_0_0(26) * occ_func_0_0(10) + occ_func_0_0(28) * occ_func_0_0(12) + occ_func_0_0(4) * occ_func_0_0(9) + occ_func_0_0(1) * occ_func_0_0(25) + occ_func_0_0(4) * occ_func_0_0(28) + occ_func_0_0(1) * occ_func_0_0(12) + occ_func_0_0(25) * occ_func_0_0(9) + occ_func_0_0(27) * occ_func_0_0(11) + occ_func_0_0(3) * occ_func_0_0(10) + occ_func_0_0(2) * occ_func_0_0(26)) / 4.;
  }

  /**** Basis functions for orbit 19****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 3.0000000 0.0000000 A  B  

1.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_19_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(8) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(7) + occ_func_0_0(8) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(8)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_19_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(8) + occ_func_0_0(26) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(5) * occ_func_0_0(16) * occ_func_0_0(0) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(7) * occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(8) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(5) + occ_func_0_0(7) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(6) + occ_func_0_0(26) * occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(8) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(14)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_19_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(27) * occ_func_0_0(8) + occ_func_0_0(26) * occ_func_0_0(17) + occ_func_0_0(5) * occ_func_0_0(16) + occ_func_0_0(28) * occ_func_0_0(20) + occ_func_0_0(7) * occ_func_0_0(13) + occ_func_0_0(25) * occ_func_0_0(6) + occ_func_0_0(28) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(20) + occ_func_0_0(25) * occ_func_0_0(13) + occ_func_0_0(8) * occ_func_0_0(17) + occ_func_0_0(27) * occ_func_0_0(16) + occ_func_0_0(26) * occ_func_0_0(5) + occ_func_0_0(7) * occ_func_0_0(18) + occ_func_0_0(27) * occ_func_0_0(6) + occ_func_0_0(26) * occ_func_0_0(15) + occ_func_0_0(27) * occ_func_0_0(18) + occ_func_0_0(26) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(15) + occ_func_0_0(28) * occ_func_0_0(19) + occ_func_0_0(8) * occ_func_0_0(14) + occ_func_0_0(25) * occ_func_0_0(5) + occ_func_0_0(28) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(19) + occ_func_0_0(25) * occ_func_0_0(14)) / 8.;
  }

  /**** Basis functions for orbit 20****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 3.0000000 0.0000000 A  B  

2.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_20_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(20) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(17) + occ_func_0_0(20) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(18)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_20_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(20) + occ_func_0_0(26) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(13) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(17) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(20) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(22) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(13) + occ_func_0_0(19) * occ_func_0_0(24) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(27) * occ_func_0_0(14) + occ_func_0_0(26) * occ_func_0_0(0) * occ_func_0_0(21) + occ_func_0_0(27) * occ_func_0_0(0) * occ_func_0_0(24) + occ_func_0_0(0) * occ_func_0_0(26) * occ_func_0_0(19) + occ_func_0_0(14) * occ_func_0_0(21) * occ_func_0_0(0) + occ_func_0_0(28) * occ_func_0_0(0) * occ_func_0_0(23) + occ_func_0_0(18) * occ_func_0_0(22) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(25) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(28) * occ_func_0_0(18) + occ_func_0_0(15) * occ_func_0_0(23) * occ_func_0_0(0) + occ_func_0_0(25) * occ_func_0_0(0) * occ_func_0_0(22)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_20_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(27) * occ_func_0_0(20) + occ_func_0_0(26) * occ_func_0_0(23) + occ_func_0_0(13) * occ_func_0_0(22) + occ_func_0_0(28) * occ_func_0_0(24) + occ_func_0_0(17) * occ_func_0_0(21) + occ_func_0_0(25) * occ_func_0_0(16) + occ_func_0_0(28) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(24) + occ_func_0_0(25) * occ_func_0_0(21) + occ_func_0_0(20) * occ_func_0_0(23) + occ_func_0_0(27) * occ_func_0_0(22) + occ_func_0_0(26) * occ_func_0_0(13) + occ_func_0_0(19) * occ_func_0_0(24) + occ_func_0_0(27) * occ_func_0_0(14) + occ_func_0_0(26) * occ_func_0_0(21) + occ_func_0_0(27) * occ_func_0_0(24) + occ_func_0_0(26) * occ_func_0_0(19) + occ_func_0_0(14) * occ_func_0_0(21) + occ_func_0_0(28) * occ_func_0_0(23) + occ_func_0_0(18) * occ_func_0_0(22) + occ_func_0_0(25) * occ_func_0_0(15) + occ_func_0_0(28) * occ_func_0_0(18) + occ_func_0_0(15) * occ_func_0_0(23) + occ_func_0_0(25) * occ_func_0_0(22)) / 8.;
  }

  /**** Basis functions for orbit 21****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

1.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_21_0() const {
    return occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(4) * occ_func_0_0(8);
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_21_0_at_0() const {
    return occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(4) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(0);
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_21_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * occ_func_0_0(3) * occ_func_0_0(4) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(3) + occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(2);
  }

  /**** Basis functions for orbit 22****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_22_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(12) * occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_22_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(9) * occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(5) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(7) * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(9)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_22_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(7) * occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(2) * occ_func_0_0(4) + occ_func_0_0(9) * occ_func_0_0(5) * occ_func_0_0(1) + occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(4) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(7) * occ_func_0_0(2) + occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(6) * occ_func_0_0(3) + occ_func_0_0(3) * occ_func_0_0(1) * occ_func_0_0(2) + occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(8) * occ_func_0_0(4) + occ_func_0_0(4) * occ_func_0_0(3) * occ_func_0_0(1) + occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(5) + occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(9)) / 4.;
  }

  /**** Basis functions for orbit 23****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

1.0000000 1.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_23_0() const {
    return occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(8) * occ_func_0_0(12);
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_23_0_at_0() const {
    return occ_func_0_0(0) * occ_func_0_0(7) * occ_func_0_0(8) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(9) * occ_func_0_0(5) * occ_func_0_0(6) * occ_func_0_0(0);
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_23_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * occ_func_0_0(7) * occ_func_0_0(8) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(11) * occ_func_0_0(8) + occ_func_0_0(5) * occ_func_0_0(10) * occ_func_0_0(7) + occ_func_0_0(9) * occ_func_0_0(5) * occ_func_0_0(6);
  }

  /**** Basis functions for orbit 24****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_24_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(4) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(8) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(7)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_24_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(4) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(8) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(4) * occ_func_0_0(8) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(15) * occ_func_0_0(5) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(7) + occ_func_0_0(3) * occ_func_0_0(6) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(2) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(1)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_24_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(4) + occ_func_0_0(2) * occ_func_0_0(17) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(3) + occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(20) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(4) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(1) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(3) + occ_func_0_0(4) * occ_func_0_0(8) * occ_func_0_0(2) + occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(5) + occ_func_0_0(4) * occ_func_0_0(19) * occ_func_0_0(7) + occ_func_0_0(3) * occ_func_0_0(6) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(2) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(1)) / 4.;
  }

  /**** Basis functions for orbit 25****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

0.0000000 2.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_25_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(11) * occ_func_0_0(4) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(12) * occ_func_0_0(3)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_25_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(11) * occ_func_0_0(4) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) + occ_func_0_0(10) * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(8) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(9) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(12) * occ_func_0_0(2) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5) + occ_func_0_0(9) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(4) * occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(11) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(10) * occ_func_0_0(1) + occ_func_0_0(4) * occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(11) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(5) + occ_func_0_0(10) * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(11) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(8) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(10) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(7) + occ_func_0_0(3) * occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(9) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(12) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6) + occ_func_0_0(9) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(14)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_25_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(11) * occ_func_0_0(4) + occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(7) + occ_func_0_0(10) * occ_func_0_0(2) * occ_func_0_0(17) + occ_func_0_0(1) * occ_func_0_0(6) * occ_func_0_0(16) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(9) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(5) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(8) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(12) * occ_func_0_0(2) + occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_0(5) + occ_func_0_0(9) * occ_func_0_0(1) * occ_func_0_0(13) + occ_func_0_0(4) * occ_func_0_0(7) * occ_func_0_0(17) + occ_func_0_0(11) * occ_func_0_0(3) * occ_func_0_0(16) + occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_0(6) + occ_func_0_0(2) * occ_func_0_0(10) * occ_func_0_0(1) + occ_func_0_0(4) * occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(3) * occ_func_0_0(11) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(3) * occ_func_0_0(5) + occ_func_0_0(10) * occ_func_0_0(2) * occ_func_0_0(15) + occ_func_0_0(11) * occ_func_0_0(3) * occ_func_0_0(18) + occ_func_0_0(3) * occ_func_0_0(2) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(10) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(5) * occ_func_0_0(15) + occ_func_0_0(12) * occ_func_0_0(4) * occ_func_0_0(19) + occ_func_0_0(4) * occ_func_0_0(1) * occ_func_0_0(7) + occ_func_0_0(3) * occ_func_0_0(6) * occ_func_0_0(14) + occ_func_0_0(1) * occ_func_0_0(9) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(12) * occ_func_0_0(3) + occ_func_0_0(2) * occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(1) * occ_func_0_0(4) * occ_func_0_0(6) + occ_func_0_0(9) * occ_func_0_0(1) * occ_func_0_0(14)) / 8.;
  }

  /**** Basis functions for orbit 26****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

1.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_26_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(3) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(4) + occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(12)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_26_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(4) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(3) + occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(6) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(2) + occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(1) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(6) * occ_func_0_0(5) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) * occ_func_0_0(11) + occ_func_0_0(7) * occ_func_0_0(4) * occ_func_0_0(10) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(6) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(15) * occ_func_0_0(1) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(4) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(8) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(11) * occ_func_0_0(0) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(0) * occ_func_0_0(10) + occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(12) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(5) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(9) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(9) * occ_func_0_0(0)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_26_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(8) + occ_func_0_0(2) * occ_func_0_0(17) * occ_func_0_0(4) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(11) + occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(10) + occ_func_0_0(4) * occ_func_0_0(20) * occ_func_0_0(3) + occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(12) + occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(6) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(9) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(12) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(3) * occ_func_0_0(9) + occ_func_0_0(1) * occ_func_0_0(13) * occ_func_0_0(2) + occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(11) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(10) + occ_func_0_0(3) * occ_func_0_0(16) * occ_func_0_0(1) + occ_func_0_0(2) * occ_func_0_0(6) * occ_func_0_0(5) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(11) + occ_func_0_0(7) * occ_func_0_0(4) * occ_func_0_0(10) + occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(6) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(1) + occ_func_0_0(3) * occ_func_0_0(18) * occ_func_0_0(4) + occ_func_0_0(2) * occ_func_0_0(8) * occ_func_0_0(7) + occ_func_0_0(6) * occ_func_0_0(1) * occ_func_0_0(11) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(10) + occ_func_0_0(8) * occ_func_0_0(3) * occ_func_0_0(12) + occ_func_0_0(4) * occ_func_0_0(19) * occ_func_0_0(2) + occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(5) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(9) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(12) + occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(5) * occ_func_0_0(2) * occ_func_0_0(9)) / 8.;
  }

  /**** Basis functions for orbit 27****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

2.0000000 -1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_27_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_27_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(6) * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(14) * occ_func_0_0(5) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(15) * occ_func_0_0(7) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(6) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(18) * occ_func_0_0(6) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(5) * occ_func_0_0(15) + occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(5) + occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(14)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_27_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(17) * occ_func_0_0(7) * occ_func_0_0(19) + occ_func_0_0(16) * occ_func_0_0(3) * occ_func_0_0(8) + occ_func_0_0(6) * occ_func_0_0(2) * occ_func_0_0(4) + occ_func_0_0(14) * occ_func_0_0(5) * occ_func_0_0(1) + occ_func_0_0(20) * occ_func_0_0(8) * occ_func_0_0(18) + occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(3) + occ_func_0_0(15) * occ_func_0_0(7) * occ_func_0_0(2) + occ_func_0_0(13) * occ_func_0_0(1) * occ_func_0_0(6) + occ_func_0_0(20) * occ_func_0_0(4) * occ_func_0_0(7) + occ_func_0_0(18) * occ_func_0_0(6) * occ_func_0_0(3) + occ_func_0_0(8) * occ_func_0_0(1) * occ_func_0_0(2) + occ_func_0_0(13) * occ_func_0_0(5) * occ_func_0_0(15) + occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(4) + occ_func_0_0(7) * occ_func_0_0(3) * occ_func_0_0(1) + occ_func_0_0(17) * occ_func_0_0(2) * occ_func_0_0(5) + occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(14)) / 4.;
  }

  /**** Basis functions for orbit 28****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_28_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(19) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(8) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_28_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(18) + occ_func_0_0(6) * occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(9) * occ_func_0_0(15) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(19) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(19) + occ_func_0_0(8) * occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(14) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(9)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_28_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(17) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(16) * occ_func_0_0(3) * occ_func_0_0(18) + occ_func_0_0(6) * occ_func_0_0(2) * occ_func_0_0(8) + occ_func_0_0(9) * occ_func_0_0(15) * occ_func_0_0(5) + occ_func_0_0(20) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(19) * occ_func_0_0(7) + occ_func_0_0(5) * occ_func_0_0(4) * occ_func_0_0(6) + occ_func_0_0(13) * occ_func_0_0(1) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(4) * occ_func_0_0(19) + occ_func_0_0(8) * occ_func_0_0(1) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(14) * occ_func_0_0(6) + occ_func_0_0(13) * occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(8) + occ_func_0_0(7) * occ_func_0_0(3) * occ_func_0_0(5) + occ_func_0_0(17) * occ_func_0_0(2) * occ_func_0_0(15) + occ_func_0_0(16) * occ_func_0_0(6) * occ_func_0_0(9)) / 4.;
  }

  /**** Basis functions for orbit 29****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_29_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(18) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(12) * occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(12) * occ_func_0_0(20) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(17) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(8)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_29_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(19) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(9) * occ_func_0_0(14) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(18) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(15) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(15) + occ_func_0_0(12) * occ_func_0_0(19) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(14) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(6) * occ_func_0_0(9) + occ_func_0_0(12) * occ_func_0_0(20) * occ_func_0_0(7) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) * occ_func_0_0(6) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(9) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(15) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(20) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(8) * occ_func_0_0(12) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(0) * occ_func_0_0(7) + occ_func_0_0(9) * occ_func_0_0(13) * occ_func_0_0(6) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(17) + occ_func_0_0(11) * occ_func_0_0(16) * occ_func_0_0(8) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(10) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(5) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(8) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(5) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(14) * occ_func_0_0(16)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_29_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(12) + occ_func_0_0(2) * occ_func_0_0(17) * occ_func_0_0(19) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(8) + occ_func_0_0(9) * occ_func_0_0(14) * occ_func_0_0(5) + occ_func_0_0(4) * occ_func_0_0(20) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(15) * occ_func_0_0(7) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(6) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(7) + occ_func_0_0(11) * occ_func_0_0(18) * occ_func_0_0(6) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(13) * occ_func_0_0(15) + occ_func_0_0(12) * occ_func_0_0(19) * occ_func_0_0(8) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(5) + occ_func_0_0(3) * occ_func_0_0(16) * occ_func_0_0(14) + occ_func_0_0(2) * occ_func_0_0(6) * occ_func_0_0(9) + occ_func_0_0(12) * occ_func_0_0(20) * occ_func_0_0(7) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(6) + occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(9) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(13) + occ_func_0_0(3) * occ_func_0_0(18) * occ_func_0_0(20) + occ_func_0_0(2) * occ_func_0_0(8) * occ_func_0_0(12) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(7) + occ_func_0_0(9) * occ_func_0_0(13) * occ_func_0_0(6) + occ_func_0_0(4) * occ_func_0_0(19) * occ_func_0_0(17) + occ_func_0_0(11) * occ_func_0_0(16) * occ_func_0_0(8) + occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(10) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(5) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(8) + occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(11) + occ_func_0_0(10) * occ_func_0_0(17) * occ_func_0_0(5) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(16)) / 8.;
  }

  /**** Basis functions for orbit 30****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

2.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_30_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(12) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(11)) / 2.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_30_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(12) * occ_func_0_0(20) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(12) + occ_func_0_0(9) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(13) * occ_func_0_0(9) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(11) + occ_func_0_0(17) * occ_func_0_0(10) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(10) * occ_func_0_0(15) * occ_func_0_0(0) * occ_func_0_0(1)) / 2.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_30_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(12) * occ_func_0_0(20) + occ_func_0_0(2) * occ_func_0_0(19) * occ_func_0_0(12) + occ_func_0_0(9) * occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(13) * occ_func_0_0(9) * occ_func_0_0(2) + occ_func_0_0(4) * occ_func_0_0(18) * occ_func_0_0(11) + occ_func_0_0(17) * occ_func_0_0(10) * occ_func_0_0(4) + occ_func_0_0(1) * occ_func_0_0(11) * occ_func_0_0(16) + occ_func_0_0(10) * occ_func_0_0(15) * occ_func_0_0(1)) / 2.;
  }

  /**** Basis functions for orbit 31****
0.0000000 0.0000000 0.0000000 A  B  

1.0000000 -2.0000000 0.0000000 A  B  

1.0000000 0.0000000 0.0000000 A  B  

2.0000000 0.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_31_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(19) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(4) * occ_func_0_0(0)) / 4.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_31_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(11) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(10) * occ_func_0_0(0) * occ_func_0_0(4) + occ_func_0_0(9) * occ_func_0_0(15) * occ_func_0_0(1) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(12) * occ_func_0_0(0) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(19) * occ_func_0_0(2) * occ_func_0_0(0) + occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(9) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(12) * occ_func_0_0(19) + occ_func_0_0(11) * occ_func_0_0(14) * occ_func_0_0(3) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(0) * occ_func_0_0(2) + occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(2) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(4) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(11) * occ_func_0_0(0) * occ_func_0_0(1) + occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(10) * occ_func_0_0(15) + occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(1) * occ_func_0_0(9)) / 4.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_31_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(17) * occ_func_0_0(4) * occ_func_0_0(12) + occ_func_0_0(16) * occ_func_0_0(11) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(10) * occ_func_0_0(4) + occ_func_0_0(9) * occ_func_0_0(15) * occ_func_0_0(1) + occ_func_0_0(20) * occ_func_0_0(3) * occ_func_0_0(11) + occ_func_0_0(2) * occ_func_0_0(12) * occ_func_0_0(3) + occ_func_0_0(10) * occ_func_0_0(19) * occ_func_0_0(2) + occ_func_0_0(13) * occ_func_0_0(9) * occ_func_0_0(14) + occ_func_0_0(20) * occ_func_0_0(12) * occ_func_0_0(19) + occ_func_0_0(11) * occ_func_0_0(14) * occ_func_0_0(3) + occ_func_0_0(3) * occ_func_0_0(9) * occ_func_0_0(2) + occ_func_0_0(13) * occ_func_0_0(2) * occ_func_0_0(10) + occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(4) + occ_func_0_0(4) * occ_func_0_0(11) * occ_func_0_0(1) + occ_func_0_0(17) * occ_func_0_0(10) * occ_func_0_0(15) + occ_func_0_0(16) * occ_func_0_0(1) * occ_func_0_0(9)) / 4.;
  }

  /**** Basis functions for orbit 32****
0.0000000 0.0000000 0.0000000 A  B  

0.0000000 1.0000000 0.0000000 A  B  

1.0000000 -1.0000000 0.0000000 A  B  

2.0000000 1.0000000 0.0000000 A  B  

  ****/
  template<typename Scalar>
  Scalar square_binary_Clexulator_default::eval_bfunc_32_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(11) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(20) * occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(19) * occ_func_0_0(12) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(12) + occ_func_0_0(18) * occ_func_0_0(11) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(20)) / 8.;
  }

  template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_eval_bfunc_32_0_at_0() const {
    return (occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(20) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(0) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(9) * occ_func_0_0(15) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(20) * occ_func_0_0(11) + occ_func_0_0(17) * occ_func_0_0(10) * occ_func_0_0(19) * occ_func_0_0(0) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(16) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(0) * occ_func_0_0(14) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(0) * occ_func_0_0(19) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(11) * occ_func_0_0(14) * occ_func_0_0(0) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(13) * occ_func_0_0(10) + occ_func_0_0(20) * occ_func_0_0(12) * occ_func_0_0(18) * occ_func_0_0(0) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(0) * occ_func_0_0(15) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(16) * occ_func_0_0(9) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(6) * occ_func_0_0(13) + occ_func_0_0(19) * occ_func_0_0(12) * occ_func_0_0(17) * occ_func_0_0(0) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(0) * occ_func_0_0(16) + occ_func_0_0(0) * occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(14) + occ_func_0_0(2) * occ_func_0_0(0) * occ_func_0_0(15) * occ_func_0_0(9) + occ_func_0_0(3) * occ_func_0_0(0) * occ_func_0_0(18) * occ_func_0_0(12) + occ_func_0_0(0) * occ_func_0_0(2) * occ_func_0_0(8) * occ_func_0_0(19) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(0) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(9) * occ_func_0_0(16) * occ_func_0_0(0) + occ_func_0_0(18) * occ_func_0_0(11) * occ_func_0_0(20) * occ_func_0_0(0) + occ_func_0_0(4) * occ_func_0_0(0) * occ_func_0_0(19) * occ_func_0_0(10) + occ_func_0_0(0) * occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(15) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(0) * occ_func_0_0(13) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(0) * occ_func_0_0(20) + occ_func_0_0(0) * occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(0) * occ_func_0_0(14) * occ_func_0_0(11) + occ_func_0_0(15) * occ_func_0_0(10) * occ_func_0_0(13) * occ_func_0_0(0)) / 8.;
  }

    template<typename Scalar>
  Scalar square_binary_Clexulator_default::site_deval_bfunc_32_0_at_0(int occ_i, int occ_f) const {
    return (m_occ_func_0_0[occ_f] - m_occ_func_0_0[occ_i]) * (occ_func_0_0(3) * occ_func_0_0(7) * occ_func_0_0(20) + occ_func_0_0(2) * occ_func_0_0(17) * occ_func_0_0(12) + occ_func_0_0(6) * occ_func_0_0(16) * occ_func_0_0(18) + occ_func_0_0(13) * occ_func_0_0(9) * occ_func_0_0(15) + occ_func_0_0(4) * occ_func_0_0(20) * occ_func_0_0(11) + occ_func_0_0(17) * occ_func_0_0(10) * occ_func_0_0(19) + occ_func_0_0(1) * occ_func_0_0(8) * occ_func_0_0(16) + occ_func_0_0(5) * occ_func_0_0(13) * occ_func_0_0(14) + occ_func_0_0(8) * occ_func_0_0(20) * occ_func_0_0(19) + occ_func_0_0(4) * occ_func_0_0(5) * occ_func_0_0(17) + occ_func_0_0(16) * occ_func_0_0(11) * occ_func_0_0(14) + occ_func_0_0(1) * occ_func_0_0(13) * occ_func_0_0(10) + occ_func_0_0(20) * occ_func_0_0(12) * occ_func_0_0(18) + occ_func_0_0(7) * occ_func_0_0(17) * occ_func_0_0(15) + occ_func_0_0(3) * occ_func_0_0(16) * occ_func_0_0(9) + occ_func_0_0(2) * occ_func_0_0(6) * occ_func_0_0(13) + occ_func_0_0(19) * occ_func_0_0(12) * occ_func_0_0(17) + occ_func_0_0(8) * occ_func_0_0(18) * occ_func_0_0(16) + occ_func_0_0(3) * occ_func_0_0(5) * occ_func_0_0(14) + occ_func_0_0(2) * occ_func_0_0(15) * occ_func_0_0(9) + occ_func_0_0(3) * occ_func_0_0(18) * occ_func_0_0(12) + occ_func_0_0(2) * occ_func_0_0(8) * occ_func_0_0(19) + occ_func_0_0(5) * occ_func_0_0(15) * occ_func_0_0(17) + occ_func_0_0(14) * occ_func_0_0(9) * occ_func_0_0(16) + occ_func_0_0(18) * occ_func_0_0(11) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(19) * occ_func_0_0(10) + occ_func_0_0(1) * occ_func_0_0(7) * occ_func_0_0(15) + occ_func_0_0(6) * occ_func_0_0(14) * occ_func_0_0(13) + occ_func_0_0(7) * occ_func_0_0(19) * occ_func_0_0(20) + occ_func_0_0(4) * occ_func_0_0(6) * occ_func_0_0(18) + occ_func_0_0(1) * occ_func_0_0(14) * occ_func_0_0(11) + occ_func_0_0(15) * occ_func_0_0(10) * occ_func_0_0(13)) / 8.;
  }

} // namespace clexulator
} // namespace CASM


extern "C" {
  /// \brief Returns a clexulator::BaseClexulator* owning a square_binary_Clexulator_default
  CASM::clexulator::BaseClexulator *make_square_binary_Clexulator_default() {
    return new CASM::clexulator::square_binary_Clexulator_default();
  }

}

