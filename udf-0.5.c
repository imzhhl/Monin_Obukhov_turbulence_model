/*周期送风边界条件*/


#include "udf.h"
#define PI 3.14
#define A 1
#define P 0.5

DEFINE_PROFILE(left_vent,thread,position)
{
	face_t f;
	begin_f_loop(f, thread)
	{
		real t = RP_Get_Real("flow-time");
		F_PROFILE(f, thread, position) = (A*0.455) *sin(2*PI*t/P)+(A*0.455);
	}
	end_f_loop(f,thread)
}

DEFINE_PROFILE(right_vent,thread,position)
{
	face_t f;
	begin_f_loop(f, thread)
	{
		real t = RP_Get_Real("flow-time");
		F_PROFILE(f, thread, position) = (A*0.455)*sin(2*PI*t/P - PI)+(A*0.455);
	}
	end_f_loop(f,thread)
}

/*源项扩散系数*/
DEFINE_DIFFUSIVITY(diff_2,c,t,i)
{
	return 1.0 * 2.88e-05 + C_MU_EFF(c,t) / 0.7;
}

DEFINE_DIFFUSIVITY(mean_age_diff, c, t, i)
{
	return C_MU_L(c,t) + C_MU_T(c, t);
}

DEFINE_SOURCE(mean_age_source, c, t, ds, eqn)
{
	ds[eqn] = 0.0;
	return C_R(c, t);
}

DEFINE_UDS_FLUX(adjoint_flux,f,t,i)
{
	if(i == 0 || i == 1)
		return F_FLUX(f, t);
	else
		return -1.0*F_FLUX(f, t);
}

