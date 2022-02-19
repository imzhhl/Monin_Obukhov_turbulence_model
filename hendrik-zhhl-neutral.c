#include "udf.h"
#include "math.h"

/* **************************************************************
** Neutral **
**************************************************************
Fluent UDFs for simulating neutral ABL flow

Control via the defined parameters
Ensure the solver is in expert mode
Use compiled UDF method
Any changes by ZHHL have been marked

Model Axis. yz = inlet/outlet plane, xz = sides, z = AGL, origin at the inlet, positive in direction of flow and AGL, fllowing along x- to x+.

C_UDMI - 3 User memory slots, 1User scalar slot
Wall Distance
Cor x
Cor y

C_UDSI
wallPhi - See description in define cell wall distance

-------------------------------------------
Owner: Hendri Breedt <u10028422@tuks.co.za>
Check: Hongliang Zhang <zhhl_email@qq.com>
Publish Date: 09/11/2017
Check Date: 01/11/2022
Version: 00 - Public release */

/* Model Constants - DTU */
#define Cmu 0.09 /*NOT 0.03!*/
#define vonKarman 0.4
#define Ce1 1.21
#define Ce2 1.92
#define sigma_k 1.0
#define sigma_e 1.3
#define sigma_theta 1.0
#define PrTurb 0.85


/* Wind speed relations */
#define z0 0.0128 /*ground roughness hight, m*/
#define Cs 0.5 /* Roughness Constant */
#define uStar 0.320 /* fraction velocity, uStar = (vonKarman*uRef)/log(zRef/z0); */
#define ablHeight 1000.0 /* Height of ABL, this is the height for fixed values of all profiles and sources */

/* General */
#define pi 3.141592
#define g -9.80665
#define R 8.3144598 /* Universal Gas Constant - Dry Air */
#define M 0.0289644 /* Molar mass of Earth's air */
#define Lb -0.0065 /* Standard temperature lapse rate */

/* Operating Conditions - Material Air */
#define presOper 101325 /* Operating Pressure Pa - Internal Solver Pressure. This is the pressure specified at 0m and for this you can use lowest mast pressure reading */
#define tempOper 288.16 /* Operating Temperature - Internal Solver Standard Tempearture. This is the temperature based from the lowest measurement height on the mast. But can be left as the standard value */
#define densOper 1.0 /* Problem density. more usual to set-1.0 ZHHL, origial is 1.0919 */  /*ZHHL*/
#define Cp 1006.43
#define beta 0.032
#define viscosity 1.7894e-05

/* Initilization */
/* Due to HAGL variations and Fluent not being able to compute cell distance before initialiazing we have to manually set the initialiaze values. These are used for z values lower than maxZInit, afterwards it returns to the inlt profile values */
#define maxZInit 0 /* Height before using init values from inlet profiles */
#define initVelocity 10.0 /* y velocity */
#define initK 2.0 /* k */
#define initEpsilon 2.0 /* epsilon */

/* ********************** Inlet Velocity ********************** */
DEFINE_PROFILE(inlet_V_Neutral, t, i)
{
	double x[ND_ND];
	double z;
	face_t f;

	begin_f_loop(f, t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0;
		if (z > ablHeight){
			z = ablHeight;
		}
		F_PROFILE(f, t, i) = (uStar/vonKarman)*log(z/z0);
	}
	end_f_loop(f, t)
}

/* ********************** Inlet k ********************** */
DEFINE_PROFILE(inlet_k_Neutral, t, i)
{
	double x[ND_ND];
	face_t f;

	begin_f_loop(f, t)
	{
		F_CENTROID(x,f,t);
		F_PROFILE(f, t, i) = pow(uStar,2.0)/sqrt(Cmu);
	}
	end_f_loop(f, t)
}

/* ********************** Inlet epsilon ********************** */
DEFINE_PROFILE(inlet_e_Neutral, t, i)
{
	double x[ND_ND];
	double z;
	face_t f;

	begin_f_loop(f, t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0;
		if (z > ablHeight){
			z = ablHeight;
		}
		F_PROFILE(f, t, i) = pow(uStar,3.0)/(vonKarman*z);
	}
	end_f_loop(f, t)
}


/* ********************** Wall Roughness ********************** */

/* Use this if you are using the ABL log law wall function */
DEFINE_PROFILE(wallRoughness,t,i)
{
	double x[ND_ND];
	face_t f;
	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		F_PROFILE(f,t,i) = z0; /* Use this if you are using the ABL log law wall function */
	}
	end_f_loop(f,t)
}

/* Modified wall roughness */
DEFINE_PROFILE(wallRoughnessModified,t,i)
{
	double x[ND_ND];
	face_t f;
	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		F_PROFILE(f,t,i) = 9.793*z0/Cs;
	}
	end_f_loop(f,t)
}


/* ************************* Initilization ************************** */

DEFINE_INIT(initNeutral,d)
{
	cell_t c;
	Thread *t;
	double x[ND_ND];
	double z;
	/* loop over all cell threads in the domain */
	thread_loop_c(t,d)
	{
		/* loop over all cells */
		begin_c_loop_all(c,t)
		{
			C_CENTROID(x,c,t);
			z = x[2] + z0;
			if (z > ablHeight){
				z = ablHeight;
			}

			if (z > maxZInit){
				C_U(c,t) = (uStar/vonKarman)*log(z/z0); /*x velocity */  /*ZHHL*/
				C_V(c,t) = 0.0; /* y velocity */  /*ZHHL*/
				C_W(c,t) = 0.0; /* z velocity */
				C_K(c,t) = pow(uStar,2.0)/sqrt(Cmu); /* k */
				C_D(c,t) = pow(uStar,3.0)/(vonKarman*z); /* epsilon */
				C_P(c,t) = 0.0; /*Pressure*/
			}
			else{
				C_U(c,t) =  initVelocity;  /*ZHHL*/
				C_V(c,t) = 0.0;  /*ZHHL*/
				C_W(c,t) = 0.0;
				C_K(c,t) = pow(uStar,2.0)/sqrt(Cmu);
				/* C_K(c,t) = initK; */
				C_D(c,t) = initEpsilon;
				C_P(c,t) = 0.0;
			}
		}
		end_c_loop_all(c,t)
	}
}


/* ************************* Wall Functions ************************** */

/* Designed around u/uStar = 1/K*log(z/z0) ref: Improved k-e model and wall function formulation for the RANS simulation of ABL flows, Parente et al
Removes the need for multiplying z0 by 9.73/Cs and can thus use roughness lengths directly from ABL modelling with first cell height = 2*z0*/

DEFINE_WALL_FUNCTIONS(ABL_logLaw, f, t, c0, t0, wf_ret, yPlus, Emod)
{
	double ustar_ground, E_prime, yPlus_prime, zp, dx_mag, wf_value;
	double mu=C_MU_L(c0,t0);
	double xf[ND_ND];
	double xc[ND_ND];
	double dx[ND_ND];

	F_CENTROID(xf, f, t);
	C_CENTROID(xc, c0,t0);

	dx[0] = xc[0] - xf[0];
	dx[1] = xc[1] - xf[1];
	dx[2] = xc[2] - xf[2];
	dx_mag = NV_MAG(dx);
	zp = dx_mag; /*ZHHL*/

	ustar_ground = pow((double)C_K(c0,t0),(double)0.5)*pow((double)Cmu, (double)0.25);
	E_prime = (mu/densOper)/(z0*ustar_ground);
	yPlus_prime = (zp+z0)*ustar_ground/(mu/densOper);

	switch (wf_ret)
	{
	case UPLUS_LAM:
		wf_value = yPlus;
		break;
	case UPLUS_TRB:
		wf_value = log(E_prime*yPlus_prime)/vonKarman;
		/*wf_value = log(Emod*yPlus)/vonKarman; Standard Fluent*/
		break;
	case DUPLUS_LAM:
		wf_value = 1.0;
		break;
	case DUPLUS_TRB:
		wf_value = 1.0/(vonKarman*yPlus_prime);
		break;
	case D2UPLUS_TRB:
		wf_value = -1.0/(vonKarman*yPlus_prime*yPlus_prime);
		break;
	default:
		printf("Wall function return value unavailable\n");
	}
	return wf_value;
}
