// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// expWArma
List expWArma(const Rcpp::List& iz_, const Rcpp::List& b_, const arma::mat& gam, const Rcpp::List& h_);
RcppExport SEXP _joineRML_expWArma(SEXP iz_SEXP, SEXP b_SEXP, SEXP gamSEXP, SEXP h_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type iz_(iz_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type b_(b_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h_(h_SEXP);
    rcpp_result_gen = Rcpp::wrap(expWArma(iz_, b_, gam, h_));
    return rcpp_result_gen;
END_RCPP
}
// gammaUpdate_approx
List gammaUpdate_approx(const Rcpp::List& b_, const Rcpp::List& z_, const Rcpp::List& w_, const Rcpp::List& pb_, const arma::vec& haz, const Rcpp::List& v_, const Rcpp::List& h_, const int& K, const int& q, const int& nev);
RcppExport SEXP _joineRML_gammaUpdate_approx(SEXP b_SEXP, SEXP z_SEXP, SEXP w_SEXP, SEXP pb_SEXP, SEXP hazSEXP, SEXP v_SEXP, SEXP h_SEXP, SEXP KSEXP, SEXP qSEXP, SEXP nevSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type b_(b_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type w_(w_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pb_(pb_SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type haz(hazSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type v_(v_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h_(h_SEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const int& >::type nev(nevSEXP);
    rcpp_result_gen = Rcpp::wrap(gammaUpdate_approx(b_, z_, w_, pb_, haz, v_, h_, K, q, nev));
    return rcpp_result_gen;
END_RCPP
}
// gammaUpdate
List gammaUpdate(const Rcpp::List& b_, const Rcpp::List& z_, const Rcpp::List& w_, const Rcpp::List& pb_, const arma::vec& haz, const Rcpp::List& v_, const Rcpp::List& h_, const int& K, const int& q, const int& nev, const arma::vec& jcount);
RcppExport SEXP _joineRML_gammaUpdate(SEXP b_SEXP, SEXP z_SEXP, SEXP w_SEXP, SEXP pb_SEXP, SEXP hazSEXP, SEXP v_SEXP, SEXP h_SEXP, SEXP KSEXP, SEXP qSEXP, SEXP nevSEXP, SEXP jcountSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type b_(b_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type z_(z_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type w_(w_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pb_(pb_SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type haz(hazSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type v_(v_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h_(h_SEXP);
    Rcpp::traits::input_parameter< const int& >::type K(KSEXP);
    Rcpp::traits::input_parameter< const int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const int& >::type nev(nevSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type jcount(jcountSEXP);
    rcpp_result_gen = Rcpp::wrap(gammaUpdate(b_, z_, w_, pb_, haz, v_, h_, K, q, nev, jcount));
    return rcpp_result_gen;
END_RCPP
}
// hazHat
arma::mat hazHat(const Rcpp::List& w_, const Rcpp::List& pb_, const arma::vec& nev);
RcppExport SEXP _joineRML_hazHat(SEXP w_SEXP, SEXP pb_SEXP, SEXP nevSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type w_(w_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pb_(pb_SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nev(nevSEXP);
    rcpp_result_gen = Rcpp::wrap(hazHat(w_, pb_, nev));
    return rcpp_result_gen;
END_RCPP
}
// lambdaUpdate
arma::mat lambdaUpdate(const Rcpp::List& b_, const Rcpp::List& imat_, const Rcpp::List& zt_, const Rcpp::List& pb_, const Rcpp::List& v_, const arma::mat& gam, const arma::vec& gam_vec, const int& q, const arma::vec& nev, const Rcpp::List& h_);
RcppExport SEXP _joineRML_lambdaUpdate(SEXP b_SEXP, SEXP imat_SEXP, SEXP zt_SEXP, SEXP pb_SEXP, SEXP v_SEXP, SEXP gamSEXP, SEXP gam_vecSEXP, SEXP qSEXP, SEXP nevSEXP, SEXP h_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type b_(b_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type imat_(imat_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type zt_(zt_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type pb_(pb_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type v_(v_SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type gam_vec(gam_vecSEXP);
    Rcpp::traits::input_parameter< const int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nev(nevSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type h_(h_SEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaUpdate(b_, imat_, zt_, pb_, v_, gam, gam_vec, q, nev, h_));
    return rcpp_result_gen;
END_RCPP
}
// mvrnormArma
arma::mat mvrnormArma(const int& n, const arma::vec& mu, const arma::mat& sigma);
RcppExport SEXP _joineRML_mvrnormArma(SEXP nSEXP, SEXP muSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(mvrnormArma(n, mu, sigma));
    return rcpp_result_gen;
END_RCPP
}
// bSim
List bSim(const int& n, const Rcpp::List& Mean_, const Rcpp::List& Sigma_);
RcppExport SEXP _joineRML_bSim(SEXP nSEXP, SEXP Mean_SEXP, SEXP Sigma_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Mean_(Mean_SEXP);
    Rcpp::traits::input_parameter< const Rcpp::List& >::type Sigma_(Sigma_SEXP);
    rcpp_result_gen = Rcpp::wrap(bSim(n, Mean_, Sigma_));
    return rcpp_result_gen;
END_RCPP
}
