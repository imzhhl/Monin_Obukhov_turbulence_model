#include "udf.h"
#include "math.h"

/* **************************************************************
** Unstable **
**************************************************************

Fluent UDFs for simulating unstable ABL flow

Control via the defined parameters
Ensure the solver is in expert mode
Use compiled UDF method

/* C_UDMI - 12 User memory slots, 1 User scalar slot
Wall Distance
Cor x
Cor y
k DTU
k Dtu Norm
epsilon Fluent
epsilon AM
epsilon AM Ce3
epsilon AM Gb
epsilon DTU
epsilon DTU - Ce3
DTU Gb

C_UDSI
wallPhi - See description in define cell wall distance

-------------------------------------------
Owner: Hendri Breedt <u10028422@tuks.co.za>
Date: 09/11/2017
Version: 00 - Public release */

/* Model Constants - DTU */
#define Cmu 0.09
#define vonKarman 0.40
#define Ce1 1.21
#define Ce2 1.92
#define sigma_k 1.0
#define sigma_e 1.3
#define sigma_theta 1.0
#define PrTurb 0.85

/* Wind speed relations */
#define z0 0.0128/*m*/
#define Cs 0.5 /* Roughness Constant */
#define uStar 0.320
#define Lin -260/* L at the inlet - L must be < 0 to use this UDF set!!! */
#define ablHeight 2000

/* General */
#define pi 3.141592
#define g -9.80665
#define R 8.3144598 /* Universal Gas Constant - Dry Air */
#define M 0.0289644 /* Molar mass of Earth's air */
#define Lb -0.0065 /* Standard temperature lapse rate*/

/* Operating Conditions - Material Air */
#define presOper 101325 /* Operating Pressure Pa - Internal Solver Pressure. This is the pressure specified at 0m and for this you can use lowest mast pressure reading */
#define tempOper 288.16 /* Operating Temperature - Internal Solver Standard Tempearture. This is the temperature based from the lowest measurement height on the mast. But can be left as the standard value */
#define densOper 1 /* Problem density */
#define Cp 1006.43
#define beta 0.0032
#define viscosity 1.7894e-05

/* Initilization */
/* Due to HAGL variations and Fluent not being able to compute cell distance before initialiazing we have to manually set the initialiaze values. These are used for z values lower than maxZInit, afterwards it returns to the inlt profile values */
#define maxZInit 0 /* Height before using init values from inlet profiles */
#define initVelocity 10.0 /* y velocity */
#define initK 2.0 /* k */
#define initEpsilon 2.0 /* epsilon */

/* ********************* Profiles ******************************
********************** Inlet Velocity ********************** */
DEFINE_PROFILE(inlet_V_Unstable,t,i)
{
	double x[ND_ND];
	double z;
	double phiM, phiE;
	face_t f;

	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0; /*ZHHL*/
		if (z > ablHeight){
			z = ablHeight;
		}
		phiE = 1.0-(z/Lin);
		phiM = pow((double)1.0-16.0*(z/Lin),(double)-0.25);
		//F_PROFILE(f, t, i) = (uStar/vonKarman)*(log(8.0*(z/z0) * (pow((double)phiM,(double)2.0))/( pow((double)phiM+1.0,(double)2.0)*(pow((double)phiM,(double)2.0)+1.0))) -pi/2.0 + 2.0*atan(1.0/phiM));
		F_PROFILE(f, t, i) = (uStar/vonKarman)*(log(8.0*(z/z0) * (pow((double)phiM,(double)4.0))/( pow((double)phiM+1.0,(double)2.0)*(pow((double)phiM,(double)2.0)+1.0))) -pi/2.0 + 2.0*atan(1.0/phiM));
		//F_PROFILE(f, t, i) = (uStar/vonKarman)*(log(z/z0)-log(0.125*    (1.0+pow((double)phiM, (double)-2.0))    *    (pow((double)1.0+1.0/phiM,(double)2.0))      )-2.0*atan(1.0/phiM)+pi/2);
	}
	end_f_loop(f,t)
}

/* ********************* Inlet k ********************** */
DEFINE_PROFILE(inlet_k_Unstable,t,i)
{
	double x[ND_ND];
	double z;
	double phiE,phiM;
	face_t f;

	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0; /*ZHHL*/
		if (z > ablHeight){
			z = ablHeight;
		}
		phiE = 1.0-(z/Lin);
		phiM = pow(1.0-16.0*(z/Lin),-0.25);
		F_PROFILE(f,t,i) = (pow((double)uStar,(double)2.0)/sqrt(Cmu))*pow((double)phiE/phiM,(double)0.5);
	}
	end_f_loop(f,t)
}

/* ********************* Inlet epsilon ********************** */
DEFINE_PROFILE(inlet_e_Unstable,t,i)
{
	double x[ND_ND];
	double z;
	double phiE, phiM;
	face_t f;

	begin_f_loop(f,t)
	{
		F_CENTROID(x,f,t);
		z = x[2] + z0; /*ZHHL*/
		if (z > ablHeight){
			z = ablHeight;
		}
		phiE = 1.0-z/Lin;
		phiM = pow(1.0-16.0*(z/Lin),-0.25);
		F_PROFILE(f,t,i) = phiE*pow(uStar,3.0)/(vonKarman*z);
	}
	end_f_loop(f,t)
}

/* ************************ k ****************************
DTU
No energy eqaution is solved with this model*/
DEFINE_SOURCE(k_source_DTU_Unstable,c,t,dS,eqn)
{
	double fUn, phiM, phiE, phiH, CkD, source, Gb, Sk, uStarLocal;
	double x[ND_ND];
	double z, L;
	C_CENTROID(x,c,t);
	z = x[2] + z0;
	L = Lin;
	if (z > ablHeight){
		z = ablHeight;
	}


	if (N_ITER > 5) {
		phiM = pow(1.0-16.0*(z/L),-0.25);
		phiE = 1.0-(z/L);
		phiH = sigma_theta*pow(1.0-16.0*(z/L),-0.5);
		uStarLocal = pow((double)C_K(c,t),(double)0.5)*pow((double)Cmu,(double)0.25)*pow((double)phiM,(double)0.25)*pow((double)phiE,(double)-0.25);/*maybe right*/
		//uStarLocal = uStar * pow((double)phiM,(double)0.25)*pow((double)phiE,(double)-0.25);
		//uStarLocal = uStar;
		fUn = 2.0-(z/L) + 8.0*(1.0-12.0*(z/L)+7.0*pow((double)z/L,(double)2.0)) - 16.0*(z/L)*(3.0-54.0*(z/L)+35.0*pow((double)z/L,(double)2.0));
		CkD = pow((double)vonKarman,(double)2)/(sigma_k*sqrt(Cmu));
		Gb = -C_MU_T(c,t)*(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2])*((z/L)/(sigma_theta))*(phiH/pow((double)phiM,(double)2)); /* DTU Formulation */
		Sk = pow((double)uStarLocal,(double)3.0)/(vonKarman*L)*((L/z)*(phiM - phiE) - (phiH)/(sigma_theta*phiM) - 0.25*CkD*pow((double)phiM,(double)6.5)*pow((double)phiE,(double)-1.5)*fUn);

		source = -densOper*Sk + Gb;
	}
	else {
		source = 0.0; /* Only run this source after 5 iterations. The gradients can cuase divergence with an illposed initilization */
		Sk = 0.0;
	}

	dS[eqn] = 0.0;
	return source;
}


/* ************************ Epsilon **********************
Epsilon is a function of the gradients and to save these the solver needs to be in expert mode
Issue: 'solve/set/expert' in the FLUENT window, and answer YES when it asks if you want to free temporary memory

Standard Fluent buoyancy treatment for epsilon
Checking advanced buoyancy treatmnent in the viscous model box adds in the formulation below
Changes in the model is made by changing Ce3 according to the AM or DTU method
Not checking the box sets Gb = 0, this term is then re added in by the sources below. Do not check the box in the viscous box! */

/* DTU */
/* Epsilon source treatment based on an anylytical expression for Ce3 */
/* No energy eqaution is solved with this model */
DEFINE_SOURCE(epsilon_source_DTU_Unstable,c,t,dS,eqn)
{
	double x[ND_ND];
	double z, L;
	double Gb, C3e, source;
	double phiM, phiH, phiE, fe;
	C_CENTROID(x,c,t);
	z = x[2] + z0; /*ZHHL*/
	L = Lin;

	if (z > ablHeight){
		z = ablHeight;
	}


	if (N_ITER > 5) {
		phiM = pow(1.0-16.0*(z/L),-0.25);
		phiE = 1.0-(z/L);
		phiH = sigma_theta*pow(1.0-16.0*(z/L),-0.5);
		fe = pow((double)phiM,(double)2.5)*(1.0-0.75*16.0*(z/L));
		/* Gb = beta*g*C_MU_T(c,t)/PrTurb*C_T_G(c,t)[2]; */ /*Standard Fluent Gb formulation, C_MU_T = Turbulent Viscosity, PrTurb = Turbulent Prandtl number, C_T_G = [partial_T/partial_xi] */
		Gb = -C_MU_T(c,t)*(C_U_G(c,t)[2]*C_U_G(c,t)[2] + C_V_G(c,t)[2]*C_V_G(c,t)[2])*((z/L)/(sigma_theta))*(phiH/pow((double)phiM,(double)2.0)); /* DTU Formulation */
		C3e = (sigma_theta/(z/L))*(phiM/phiH)*(Ce1*phiM-Ce2*phiE+(Ce2-Ce1)*pow((double)phiE,(double)-0.5)*fe); /* DTU C3e formulation */

		source = Ce1*C_D(c,t)/C_K(c,t)*C3e*Gb; /*C_D = epsilon, C_K = k */
	}
	else {
		source = 0.0; /* Only run this source after 5 iterations. The gradients can cuase divergence with an illposed initilization */
		Gb = 0.0;
		C3e = 0.0;

	}
	dS[eqn] = 0.0;
	return source;
}

/* ************************* Initilization ************************** */

DEFINE_INIT(initUnstable,d)
{
	cell_t c;
	Thread *t;
	double x[ND_ND];
	double phiM, phiE, phiH, pressure, potenTemp, z, zAMSL, L;
	/* loop over all cell threads in the domain */
	thread_loop_c(t,d)
	{
		/* loop over all cells */
		begin_c_loop_all(c,t)
		{
			C_CENTROID(x,c,t);
			z = x[2] + z0;
			L = Lin;
			if (z > ablHeight){
				z = ablHeight;
			}

			if (z > maxZInit){
				phiM = pow(1.0-16.0*(z/L),-0.25);
				phiE = 1.0-(z/L);
				phiH = sigma_theta*pow(1.0-16.0*(z/L),-0.5);
				C_U(c,t) = (uStar/vonKarman)*(log(8.0*(z/z0) * (pow((double)phiM,(double)4.0))/( pow((double)phiM+1.0,(double)2.0)*(pow((double)phiM,(double)2.0)+1.0))) -pi/2.0 + 2.0*atan(1.0/phiM)); /*x velocity */ /*ZHHL*/
				C_V(c,t) =  0.0; /* y velocity */
				C_W(c,t) = 0.0; /* z velocity */
				/* C_T(c,t) = potenTemp/(pow(presOper/pressure,0.286)); /* Temperature */
				C_K(c,t) = (pow((double)uStar,(double)2.0)/sqrt(Cmu))*pow((double)phiE/phiM,(double)0.5); /* k */
				C_D(c,t) = phiE*pow(uStar,3.0)/(vonKarman*z); /* epsilon */
				C_P(c,t) = 0.0; /*Pressure*/
			}
			else{
				C_U(c,t) = 0.0;
				C_V(c,t) = initVelocity;
				C_W(c,t) = 0.0;
				C_K(c,t) = initK;
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
	zp = dx_mag/2.0; /*ZHHL*/

	ustar_ground = pow((double)C_K(c0,t0),(double)0.5)*pow((double)Cmu,(double) 0.25);
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

/* ********************** Walls *********************** 189
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