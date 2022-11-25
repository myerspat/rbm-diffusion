#include "rbm/rbm.hpp"
#include "unit_test_framework.hpp"


TEST(ContructF_t_test)
{
 xt::xarray<double> F = xt::xarray<double>::from_shape(
    {3, 3});

    F(0,0)=2.2;
    F(1,0)=0.0;
    F(2,0)=0.0;
    F(0,1)=0.0;
    F(1,1)=1.3;
    F(2,1)=0.0;
    F(0,2)=0.0;
    F(1,2)=0.0;
    F(2,2)=3.5;  

 xt::xarray<double> Flux = xt::xarray<double>::from_shape(
    {2, 3});
    
    Flux(0,0)=0.2;
    Flux(1,0)=0.94;
    Flux(0,1)=0.03;
    Flux(1,1)=0.3;
    Flux(0,2)=0.230;
    Flux(1,2)=0.540;

 // xt::xarray<double> F_t =  rbm::RBM::constructF_t(F, Flux);
//  rbm::RBM::constructF_t(F, Flux);
   rbm::PerturbAbsorption object;
   xt::xarray<double> F_t = object.constructF_t(F, Flux);

   ASSERT_ALMOST_EQUAL(F_t(0, 0), 1.2398, 0.001);

//    for (size_t i : I) {
//     for (size_t j : J) {
//       ASSERT_ALMOST_EQUAL(result.at(i, j), correct_mat[i][j], 0.000001);
//     }
//   }
    

}



TEST_MAIN();
