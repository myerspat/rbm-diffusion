#include "rbm/rbm.hpp"
#include "unit_test_framework.hpp"

TEST(ContructM_t_test)
{
 xt::xarray<double> M = xt::xarray<double>::from_shape(
    {3, 3});

    M(0,0)=2.2;
    M(1,0)=0.0;
    M(2,0)=0.0;
    M(0,1)=0.0;
    M(1,1)=1.3;
    M(2,1)=0.0;
    M(0,2)=0.0;
    M(1,2)=0.0;
    M(2,2)=3.5;  

 xt::xarray<double> Flux = xt::xarray<double>::from_shape(
    {3, 2});
    
    Flux(0,0)=0.2;
    Flux(1,0)=0.94;
    Flux(2,0)=0.03;
    Flux(0,1)=0.3;
    Flux(1,1)=0.230;
    Flux(2,1)=0.540;


   rbm::PerturbAbsorption object;
   xt::xarray<double> M_t = object.constructM_t(M, Flux);

   ASSERT_ALMOST_EQUAL(M_t(0, 0), 1.2398, 0.001);


}

TEST_MAIN();
