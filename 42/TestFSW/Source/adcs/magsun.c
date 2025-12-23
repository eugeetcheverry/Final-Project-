#include "42.h"
#include <math.h>
#include "abMATH.h"

// Constants
const double pi=3.14159265e0;
const double DEGTORAD = pi/180.;
const double RADTODEG = 180./pi;
const double RADSECTORPM = 30./pi;
const double RPMTORADSEC = pi/30.;
const double EARTHGRAVCO = 3.986004418e14; // in METER^3/SECOND^2 (es mu=GM)
const double RE = 6.378135e06;             // in METERS
const double RE2 = RE*RE;
const double EARTHOMEGA = 7.2921151467e-5; // in RADIANS/SECOND
const double EARTH_J2 =  1.082616e-3;
const double EARTH_J3 = -2.538810e-6;
const double EARTH_J4 = -1.655970e-6;

double earthradius = RE;

/* ********************************************************************
   ADCS based on TAM, CSS, Gyros, MTQ and one wheel as optional
   Detumbling: 
      TAM, Gyros and MTQ
   Sun Pointing: 
      iwheel<0 => TAM, 6xCSS, 3xGyros, and 3xMTQ
      iwheel>=0 && iwheel<3 => TAM, 6xCSS, 3Gyros, 3xMTQ and 1xWheel
   ********************************************************************/
int adcsMagSun(struct AcType *AC)
{

   struct AcCfsCtrlType *C;
   C = &AC->CfsCtrl;

   int retval = 0.;
   static int secondspointed = 0, secondssun = 0, spinseconds = 0, timepostdetumbling = 0;
   double qe[4]={0.,0.,0.,1.};

   timepostdetumbling += AC->DT;

   // Matriz Gamma promedio móvil
   static double GAV[3][3] = {{0,0,0}, {0,0,0}, {0,0,0}};

   double mbvb = MAGV(AC->bvb);        // Norm of B
   static double kq = 0., kw = 10.;    // Initial Gains
   int coriolis = 1;                   // == 1 => coriolis compensation
   int iwearth = 1;                    // >= 0 => M3, axis of rotation, x=0 y=1 z=2
   static double earthvel = 0.;        // != 0 => M3,                   earthvel * wo
   double kw0 = 20.0;                  // Derivative Gain with Eclipse
   double kw1 = 5.0;                   // Derivative Gain without Eclipse
   double kq0 = 0.075;                 // Proportional Gain with Eclipse
   double kq1 = 0.0075;                // Proportional Gain without Eclipse
   double eps = 0.001;                 // Averaging Parameter
   double thnw = 0.001;                // End of Detumbling Velocity Condition [rad/sec] 
   double nm=0.;                       // Auxiliary Variable
   double threclipse = 0.5;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
   double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);   // Norm of angular velocity in b
   int iwheel = 2; // x=0 to z=2 uses one wheel on this axis, otherwise does not use it (must be consistent with SC_xx.txt)
   if (AC->Nwhl==0) iwheel = -2;
   else if (AC->Whl[0].Axis[0]==1) iwheel = 0;
   else if (AC->Whl[0].Axis[1]==1) iwheel = 1;
   else if (AC->Whl[0].Axis[2]==1) iwheel = 2;
   else iwheel = -1;
   if (mbvb > 0.) {
      double m[3], u[3];
      for (int i = 0; i < 3; i++) {
         /* Angular Velocity Feedback */
         u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;
         if (i==iwearth && kq!=0.) {
            //if (earthvel==0) spinseconds=0; 
            //else if (spinseconds<100*60) spinseconds++;
            u[i] = u[i] + kw * eps * J[i] * 0.001082 * earthvel / mbvb / mbvb;
         }
         /* Sun Pointing Feedback */
         double sumillum = 0.;
         for (int ss=0; ss<AC->Ncss; ss++) sumillum += AC->CSS[ss].Illum;
         if (sumillum>threclipse) {
            /* Sun with AC->svb */
            secondssun++;
            double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
            double qc[4] = {-svb[2], 0., svb[0], svb[1]};
            
            static int secondschange = 0;
            static double signq = 1.;
            if (svb[1]<0. && secondschange>100*60) {
               signq = -1.;
               secondschange=0;
            } if (svb[1]>0. && secondschange>100*60) {
               signq = 1.;
               secondschange=0;
            } else {
               secondschange++;
            }
            qc[0] = signq*qc[0];
            qc[1] = signq*qc[1];
            qc[2] = signq*qc[2];
            
            if (nw < thnw && kq==0.) {kq = kq0; kw = kw0;}
            if (kq!=0.) {
               if (secondssun>2*100*60) {
                  kq = kq1;
                  kw = kw1;
               } else {
                  kq = kq0;
                  kw = kw0;
               }
            }
            if (nw < thnw && fabs(svb[1])>0.98) {
               secondspointed++;
            } else if (fabs(svb[1])<0.25) {
               secondspointed=0.;
            }
            u[i] = u[i] - kq * eps * eps * (1/J[i]) * qc[i] / mbvb / mbvb;
            if (i==iwheel && kq!=0.) {
               AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - kq * eps * eps * (1/J[i]) * qc[i];
               if (i==iwearth) AC->Tcmd[i] += kw * eps * J[i] * 0.001082 * earthvel;
               AC->Tcmd[i] = 10.0*AC->Tcmd[i];
               u[i] = 0.; // Descargar momento de ruedas
            }
            qe[0]=qc[0];qe[1]=qc[1];qe[2]=qc[2];qe[3]=qc[3];
         } else {
            if(kq!=0.) {
               if (i==iwearth) u[i] = u[i] + kw * eps * J[i] * 0.001082 * earthvel / mbvb / mbvb; 
               u[i] = 0.5*u[i];
               if (i==iwheel) {
                  AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i];
                  if (i==iwearth) AC->Tcmd[i] += kw * eps * J[i] * 0.001082 * earthvel;
                  AC->Tcmd[i] = 10.0*AC->Tcmd[i];
                  u[i] = 0.;
               }
            }
            secondspointed=0.;
            secondssun=0;
         }
      }
      /* m = b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      /* Direct Control Allocation */
      int dnm = 0;
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[0])*(1-dnm)+nm*dnm>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/(nm*dnm+(1-dnm)*fabs(m[0]));
         m[1]= m[1]*AC->MTB[0].Mmax/(nm*dnm+(1-dnm)*fabs(m[0]));
         m[2]= m[2]*AC->MTB[0].Mmax/(nm*dnm+(1-dnm)*fabs(m[0]));
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[1])*(1-dnm)+nm*dnm>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/(nm*dnm+(1-dnm)*fabs(m[1]));
         m[1]= m[1]*AC->MTB[1].Mmax/(nm*dnm+(1-dnm)*fabs(m[1]));
         m[2]= m[2]*AC->MTB[1].Mmax/(nm*dnm+(1-dnm)*fabs(m[1]));
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[2])*(1-dnm)+nm*dnm>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/(nm*dnm+(1-dnm)*fabs(m[2]));
         m[1]= m[1]*AC->MTB[2].Mmax/(nm*dnm+(1-dnm)*fabs(m[2]));
         m[2]= m[2]*AC->MTB[2].Mmax/(nm*dnm+(1-dnm)*fabs(m[2]));
      }
      /* .. Momentum Management */
      C->Kunl = 0.5E5;
      if (iwheel<2) C->Kunl = 2.5E5;
      double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
      if (fabs(m[0])<AC->MTB[0].Mmax/2. && fabs(m[1])<AC->MTB[1].Mmax/2. && fabs(m[2])<AC->MTB[2].Mmax/2.) {
         if (iwheel>=0 && iwheel<3) {
            um[iwheel] -= C->Kunl*AC->Whl[0].H;
            mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
            mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
            mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
            tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
            tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
            tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
            AC->Tcmd[iwheel] -= tm[iwheel];
         }

         double vii[3] = {0.,0.,0.}, MM[3][3], M[3];
         static int numgamma = 0;
         numgamma ++;
         for (int ii=0;ii<3;ii++) {
            vii[ii]     = 1;
            M[0]        = - ( vii[1]*AC->bvb[2] - vii[2]*AC->bvb[1] );
            M[1]        = - ( vii[2]*AC->bvb[0] - vii[0]*AC->bvb[2] );
            M[2]        = - ( vii[0]*AC->bvb[1] - vii[1]*AC->bvb[0] );
            MM[0][ii]   =   ( M[1]*AC->bvb[2] - M[2]*AC->bvb[1] ) / mbvb / mbvb;
            MM[1][ii]   =   ( M[2]*AC->bvb[0] - M[0]*AC->bvb[2] ) / mbvb / mbvb;
            MM[2][ii]   =   ( M[0]*AC->bvb[1] - M[1]*AC->bvb[0] ) / mbvb / mbvb;
            GAV[0][ii]  = 0.0001 * MM[0][ii] + 0.9999 * GAV[0][ii];
            GAV[1][ii]  = 0.0001 * MM[1][ii] + 0.9999 * GAV[1][ii];
            GAV[2][ii]  = 0.0001 * MM[2][ii] + 0.9999 * GAV[2][ii];
            vii[ii]     = 0;
         }
         
         if (coriolis==1 && kq!=0. && iwheel>0) {
            double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2]-AC->wbn[1]*AC->Whl[iwheel].H, (J[2]-J[0])*AC->wbn[2]*AC->wbn[0]+AC->wbn[0]*AC->Whl[iwheel].H, (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
            //double cor[3] = { -AC->wbn[1]*AC->Whl[iwheel].H , AC->wbn[0]*AC->Whl[iwheel].H, 0.};
            AC->Tcmd[iwheel] -= cor[iwheel];
            cor[iwheel] = 0;
            mm[0] += (cor[1]*AC->bvb[2]-cor[2]*AC->bvb[1])/0.5;
            mm[1] += (cor[2]*AC->bvb[0]-cor[0]*AC->bvb[2])/0.9;
            mm[2] += (cor[0]*AC->bvb[1]-cor[1]*AC->bvb[0])/0.5;
         }
      }
      for (int i = 0; i < 3; i++) {
         AC->Mcmd[i] = m[i] + mm[i];
      }

      /*
      if (secondssun<2500*2 && secondspointed>100*2) {
         retval = 1;
         printf("\n RETVAL=1 ss = %d, sp = %d ", secondssun, secondspointed);
      }
      */
      
      // Adaptation of drag model
      if (iwheel>=0 && kq!=0.) {
         double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
         static double chat = 0.00023/0.00177/0.00177;
         double pchat = 1E8;
         double Gamma[3] = {0.55,0.9,0.55}; 
         chat = chat + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[iwheel]*AC->wbn[iwheel]/Gamma[iwheel];
         if (chat>1.5*chat0) chat=1.15*chat0;
         if (chat<0.5*chat0) chat=0.85*chat0;
         AC->Tcmd[iwheel] -= chat*AC->Whl[0].H*fabs(AC->Whl[0].H);
         // Saturation of wheel
         double satw = 0.1;
         if ( AC->Tcmd[iwheel] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] =  AC->Whl[0].Tmax*satw;
         if ( AC->Tcmd[iwheel] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] = -AC->Whl[0].Tmax*satw;
      }

   }

   // Mode 3 Test
   // To be sent by command, this is just an example for dawn/dusk
   if ( secondssun > 30000*2 || timepostdetumbling > 10000) {
      earthvel = -1.0;
   }
   retval = 0.;

   static int first=1;
   FILE *FilePtr;
   if (first) {
      first=0;
      FilePtr = fopen("TSAT/mission.m", "w");
      fprintf(FilePtr, "vm=[%f %f %f %f %f %f %f %f %f %f %f %f %f %f];\n", qe[0], qe[1], qe[2], qe[3], AC->Tcmd[iwheel], GAV[0][0], GAV[0][1], GAV[0][2], GAV[1][0], GAV[1][1], GAV[1][2], GAV[2][0], GAV[2][1], GAV[2][2]);
      fclose(FilePtr);
   } else {
      FilePtr = fopen("TSAT/mission.m", "a");
      fprintf(FilePtr, "vm=[vm;%f %f %f %f %f %f %f %f %f %f %f %f %f %f];\n", qe[0], qe[1], qe[2], qe[3], AC->Tcmd[iwheel], GAV[0][0], GAV[0][1], GAV[0][2], GAV[1][0], GAV[1][1], GAV[1][2], GAV[2][0], GAV[2][1], GAV[2][2]);
      fclose(FilePtr);
   }

   return retval;

}


int adcsRwTriadTLE(struct AcType *AC)
{

   int sunwheel = 0;

   int retval = 0.;
   static int secondsmis = 0, secondsnmis = 0;

   secondsnmis++;
   if (secondsnmis>21000) {
      secondsmis++;
   }
   if (secondsmis>2500) {
      secondsmis = 0;
      secondsnmis = 0;
   }

   struct AcCfsCtrlType *C;
   double L1[3],L2[3],L3[3];
   double HxB[3];
   double AngErr;
   long i,j;

   C = &AC->CfsCtrl;

   if (C->Init) {
      C->Init = 0;
      for(i=0;i<3;i++) FindPDGains(AC->MOI[i][i],0.1,0.7,&C->Kr[i],&C->Kp[i]);
      C->Kunl = 1.0E6;
   }

/* .. Commanded Attitude */
   double qln2[4]={0.,0.,0.,1.}, qbr2[4]={0.,0.,0.,1.};
   if (AC->GPS[0].Valid) {
      CopyUnitV(AC->PosN,L3);
      VxV(AC->PosN,AC->VelN,L2);
      UNITV(L2);
      UNITV(L3);
      //for(i=0;i<3;i++) {
      //   L2[i] = -L2[i];
      //   L3[i] = -L3[i];
      //}
      VxV(L2,L3,L1);
      UNITV(L1);
      double phi=0.*3.1416/180.;
      if (secondsmis==0) phi=0.*3.1416/180.;
      for(i=0;i<3;i++) {
         AC->CLN[0][i] = L1[i]*cos(phi)-L3[i]*sin(phi);
         AC->CLN[1][i] = L2[i];
         AC->CLN[2][i] = L3[i]*cos(phi)+L1[i]*sin(phi);
      }
      C2Q(AC->CLN,AC->qln);
      double phi2=15.*3.1416/180.;
      if (secondsmis==0) phi2=0.*3.1416/180.;
      for(i=0;i<3;i++) {
         AC->CLN[0][i] = L1[i]*cos(phi2)-L3[i]*sin(phi2);
         AC->CLN[1][i] = L2[i];
         AC->CLN[2][i] = L3[i]*cos(phi2)+L1[i]*sin(phi2);
      }
      C2Q(AC->CLN,qln2);
      AC->wln[1] = MAGV(AC->VelN)/MAGV(AC->PosN);
   } else {
      for(i=0;i<3;i++) {
         for(j=0;j<3;j++) {
            AC->CLN[i][j] = 0.0;
         }
         AC->CLN[i][i] = 1.0;
         AC->qln[i] = 0.0;
         AC->wln[i] = 0.0;
      }
      AC->qln[3] = 1.0;
   }
            
/* .. Attitude Control Reference */
   QxQT(AC->qbn,AC->qln,AC->qbr);
   RECTIFYQ(AC->qbr);

   QxQT(AC->qbn,qln2,qbr2);
   RECTIFYQ(qbr2);

   static double sqc[4]={0.,0.,0.,0.};

   double mbvb = MAGV(AC->bvb);        // Norm of B
   static double kq = 0., kw = 10.;    // Initial Gains
   int coriolis = 0;                   // == 1 => coriolis compensation
   int iwearth = 1;                    // >= 0 => axis of rotation, x=0 y=1 z=2
   static double earthvel = 1.0;       // != 0 => earthvel * wo
   static double earthang = 0.0;       // != 0 => shift
   double kw0 = 20.0;                  // Derivative Gain with Eclipse
   double kw1 = 5.0;                   // Derivative Gain without Eclipse
   double kq0 = 0.075;                 // Proportional Gain with Eclipse
   double kq1 = 0.0075;                // Proportional Gain without Eclipse
   double eps = 0.001;                 // Averaging Parameter
   if (secondsmis>0) eps = 0.002;
   double thnw = 0.001;                // End of Detumbling Velocity Condition [rad/sec] 
   double nm=0.;                       // Auxiliary Variable
   double threclipse = 0.5;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
   double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);   // Norm of angular velocity in b
   int iwheel = 2; // x=0 to z=2 uses one wheel on this axis, otherwise does not use it (must be consistent with SC_xx.txt)
   if (AC->Nwhl==0) iwheel = -2;
   else if (AC->Whl[0].Axis[0]==1) iwheel = 0;
   else if (AC->Whl[0].Axis[1]==1) iwheel = 1;
   else if (AC->Whl[0].Axis[2]==1) iwheel = 2;
   else iwheel = -1;
   if (mbvb > 0.) {
      double m[3], u[3];
      for (int i = 0; i < 3; i++) {
         /* Angular Velocity Feedback */
         u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;
         if (i==iwearth && kq!=0.) u[i] = u[i] + kw * eps * J[i] * 0.001082 * earthvel / mbvb / mbvb;
         if (i==iwearth && kq!=0.) u[i] = u[i] + kq * eps * eps / J[i] * earthang / mbvb / mbvb;
         double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
         //double qc[4] = {AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3]};
         double qc[4] = {qbr2[0], qbr2[1], qbr2[2], qbr2[3]};
         double sqc3=1.;if (qc[3]<0.) sqc3=-1.;
         qc[0]=qc[0]*sqc3;
         qc[1]=qc[1]*sqc3;
         qc[2]=qc[2]*sqc3;
         qc[3]=qc[3]*sqc3;
         kq = kq0;
         kw = kw0;
         double ki = kq*0.;
         sqc[i] = 0.001*qc[i] + 0.999*sqc[i];
         u[i] = u[i] - kq * eps * eps * (1/J[i]) * qc[i] / mbvb / mbvb;
         u[i] = u[i] - ki * eps * eps * (1/J[i]) * sqc[i] / mbvb / mbvb;
         if (i==iwheel && kq!=0.) {
            AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - kq * eps * eps * (1/J[i]) * (qc[i]+sqc[i]*ki/kq);
            if (i==iwearth) AC->Tcmd[i] += kw * eps * J[i] * AC->wln[1] * earthvel;
            AC->Tcmd[i] = 10.0*AC->Tcmd[i];
            //u[i] = 0.; // Descargar momento de ruedas
         }
      }
      /* m = b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      /* Direct Control Allocation */
      if (fabs(m[0])>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/fabs(m[0]);
         m[1]= m[1]*AC->MTB[0].Mmax/fabs(m[0]);
         m[2]= m[2]*AC->MTB[0].Mmax/fabs(m[0]);
      }
      if (fabs(m[1])>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/fabs(m[1]);
         m[1]= m[1]*AC->MTB[1].Mmax/fabs(m[1]);
         m[2]= m[2]*AC->MTB[1].Mmax/fabs(m[1]);
      }
      if (fabs(m[2])>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/fabs(m[2]);
         m[1]= m[1]*AC->MTB[2].Mmax/fabs(m[2]);
         m[2]= m[2]*AC->MTB[2].Mmax/fabs(m[2]);
      }
      /* .. Momentum Management */
      C->Kunl = 0.5E5;
      if (iwheel<2) C->Kunl = 2.5E5;
      double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
      if (fabs(m[0])<AC->MTB[0].Mmax/2. && fabs(m[1])<AC->MTB[1].Mmax/2. && fabs(m[2])<AC->MTB[2].Mmax/2.) {
         if (iwheel>=0 && iwheel<3) {
            um[iwheel] -= C->Kunl*AC->Whl[0].H;
            mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
            mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
            mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
            tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
            tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
            tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
            AC->Tcmd[iwheel] -= tm[iwheel];
         }
         if (coriolis==1 && kq!=0.) {
            double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2], (J[2]-J[0])*AC->wbn[2]*AC->wbn[0], (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
            AC->Tcmd[iwheel] -= cor[iwheel];
            cor[iwheel]=0;
            mm[0] += (cor[1]*AC->bvb[2]-cor[2]*AC->bvb[1]);
            mm[1] += (cor[2]*AC->bvb[0]-cor[0]*AC->bvb[2]);
            mm[2] += (cor[0]*AC->bvb[1]-cor[1]*AC->bvb[0]);
         }
      }
      for (int i = 0; i < 3; i++) {
         AC->Mcmd[i] = m[i] + mm[i];
      }
      
      // Adaptation of drag model
      if (iwheel>=0 && kq!=0.) {
         double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
         static double chat = 0.00023/0.00177/0.00177;
         double pchat = 1E8;
         chat = chat + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[iwheel]*AC->wbn[iwheel];
         if (chat>1.15*chat0) chat=1.15*chat0;
         if (chat<0.85*chat0) chat=0.85*chat0;
         AC->Tcmd[iwheel] -= chat*AC->Whl[0].H*fabs(AC->Whl[0].H);
         // Saturation of wheel
         double satw = 0.1;
         if ( AC->Tcmd[iwheel] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] =  AC->Whl[0].Tmax*satw;
         if ( AC->Tcmd[iwheel] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[iwheel] = -AC->Whl[0].Tmax*satw;
      }

   }

   static int first=1;
   FILE *FilePtr;
   if (first) {
            first=0;
            FilePtr = fopen("TSAT/mission.m", "w");
            fprintf(FilePtr, "vm=[%f %f %f %f %f];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[iwheel]);
            fclose(FilePtr);
   } else {
            FilePtr = fopen("TSAT/mission.m", "a");
            fprintf(FilePtr, "vm=[vm;%f %f %f %f %f];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[iwheel]);
            fclose(FilePtr);
   }

   return retval;

}


int adcsRwStkGPS(struct AcType *AC)
{

   int retval = 0.;
   static int secondspointed = 0, secondssun = 0;

   /**
    * Detumbling y apuntamiento al Sol con CSS, TAM, gyros y bobinas
    */
   double mbvb = MAGV(AC->bvb);
   static double kq = 0.;
   double kw = 10., kw0 = 1.0, kq0 = 0.05, eps=0.001, nm=0., thnw = 0.001;
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};
   double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);
   if (mbvb > 0.) {
      double m[3],u[3];
      for (int i = 0; i < 3; i++) {
         /* Damping con gyros ideales */
         u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;
         /* Sun Pointing */
         double sumillum = 0.;
         for (int ss=0; ss<AC->Ncss; ss++) sumillum += AC->CSS[ss].Illum;
         if (sumillum>.1) {
            if (nw < thnw && kq==0.) {kq= kq0;kw=kw0;}
            /* Sol con AC->svb */
            secondssun++;
            double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
            double nsvb = sqrt(svb[0]*svb[0]+svb[1]*svb[1]+svb[2]*svb[2]);
            double qc[4] = {0., svb[2]/nsvb, -svb[1]/nsvb, 0.};
            if (svb[0]<0. && svb[0]>-1.) {
               qc[0] = 0.;
               qc[1] = 0.1*qc[1]/sqrt(qc[1]*qc[1]+qc[2]*qc[2]);
               qc[2] = 0.1*qc[2]/sqrt(qc[1]*qc[1]+qc[2]*qc[2]);
               secondspointed=0.;
            } else if (svb[0]==-1.) {
               qc[0] = 0.;
               qc[1] = 1.;
               qc[2] = 1.;
               secondspointed=0.;
            } else if (nw < thnw && svb[0]>0.98) {
               secondspointed++;
            } else {
               secondspointed=0.;
            }
            u[i] = u[i] - kq * eps * eps * (1/J[i]) * qc[i] / mbvb / mbvb;          
         } else {
            u[i] = 0.5*u[i];
            secondspointed=0.;
            secondssun=0;
         }
      }
      /* m = b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      /* Direct Control Allocation */
      if (fabs(m[0])>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/fabs(m[0]);
         m[1]= m[1]*AC->MTB[0].Mmax/fabs(m[0]);
         m[2]= m[2]*AC->MTB[0].Mmax/fabs(m[0]);
      }
      if (fabs(m[1])>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/fabs(m[1]);
         m[1]= m[1]*AC->MTB[1].Mmax/fabs(m[1]);
         m[2]= m[2]*AC->MTB[1].Mmax/fabs(m[1]);
      }
      if (fabs(m[2])>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/fabs(m[2]);
         m[1]= m[1]*AC->MTB[2].Mmax/fabs(m[2]);
         m[2]= m[2]*AC->MTB[2].Mmax/fabs(m[2]);
      }
      for (int i = 0; i < 3; i++) {
        AC->MTB[i].Mcmd = m[i];
      }
      if (secondssun<3000 && secondspointed>500) {
         retval = 1;
      }
   }

   return retval;

}

/* ********************************************************************
   FOCUS
/* ********************************************************************

/* ********************************************************************
   ADCS based on TAM, CSS, Gyros, MTQ and a Wheel pyramid
   Detumbling: 
      TAM, Gyros and MTQ
   Sun Pointing: 
      All sensors and actuators
   ********************************************************************/
int adcsMagSunf(struct AcType *AC)
{

   struct AcCfsCtrlType *C;
   C = &AC->CfsCtrl;

   int retval = 0.;
   static int secondspointed = 0, secondssun = 0;

   double mbvb = MAGV(AC->bvb);        // Norm of B

   int coriolis = 0;                   // == 1 => coriolis compensation
   static double kq = 0., kw = 10.;    // Initial Gains
   double kw0 = 5.0;                   // Derivative Gain with Eclipse
   double kw1 = 5.0;                   // Derivative Gain without Eclipse
   double kq0 = 0.075*60;              // Proportional Gain with Eclipse
   double kq1 = 0.0075*60;             // Proportional Gain without Eclipse
   double eps = 0.0015;                // Averaging Parameter
   double thnw = 0.005;             // End of Detumbling Velocity Condition [rad/sec] 
   double nm = 0.;                     // Auxiliary Variable
   double threclipse = 0.5;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
   double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);   // Norm of angular velocity in b

   if (mbvb > 0.) {
      double m[3], u[3];
      for (int i = 0; i < 3; i++) {

         /* Angular Velocity Feedback */
         u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;
       
         /* Sun Pointing Feedback */
         double sumillum = 0.;
         for (int ss=0; ss<AC->Ncss; ss++) sumillum += AC->CSS[ss].Illum;
         if (sumillum>threclipse) {
            static double sgq0=1.;
            /* Sun with AC->svb */
            double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
            double qc[4] = {-svb[2], 0., svb[0], svb[1]};
            if (i==0) {
               secondssun++;
               qc[3] = cos(acos(svb[1])/2.);
               if (svb[1]<1. && svb[1]>-1.) {
                  qc[0] *= sin(acos(svb[1])/2.)/sqrt(svb[0]*svb[0]+svb[2]*svb[2]);
                  qc[2] *= sin(acos(svb[1])/2.)/sqrt(svb[0]*svb[0]+svb[2]*svb[2]);
               } else if (svb[1]==-1.) qc[0]=1.;
               if (qc[3]<0.) sgq0=-1.;
               if (nw < thnw && kq==0.) {
                  kq = kq0; 
                  kw = kw0;
               }
               if (nw < thnw && fabs(svb[1])>0.98) {
                  secondspointed++;
               } else if (fabs(svb[1])<0.25) {
                  secondspointed=0.;
               }
            }
            if (kq!=0.) {
               AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - sgq0 * kq * eps * eps * (1/J[i]) * qc[i];
               u[i] = 0.; // Descargar momento de ruedas
               AC->Tcmd[i] = 0.1*AC->Tcmd[i]; // 10
               FindPDGains(AC->MOI[i][i],0.035,1.0, &C->Kr[i],&C->Kp[i]);
               AC->Tcmd[i] = -C->Kr[i]*AC->wbn[i]-C->Kp[i]*(2.0*qc[i]);
               double rho = 1;
               while (fabs(AC->Tcmd[i]) > AC->Whl[i].Tmax && rho>0.05) {
                  rho = rho*0.97;
                  AC->Tcmd[i] =  -rho*C->Kr[i]*AC->wbn[i] - C->Kp[i]*(2.0*(qc[i]));
               }
            } else if (nw < thnw) {
               AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i];
               AC->Tcmd[i] = 0.2*AC->Tcmd[i];
               u[i] = 0.; // Descargar momento de ruedas
            }
         } else {
            if(kq!=0.) {
               AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i];
               AC->Tcmd[i] = 0.1*AC->Tcmd[i]; // 10
               u[i] = 0.;
            }
            secondspointed=0.;
            secondssun=0;
         }

         printf("\n sumillum = %f u[%d] = %f    AC->Tcmd[%d] = %f ", sumillum, i, u[i], i, AC->Tcmd[i]);
      
      }
      /* m = b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      m[0]*=10.;m[1]*=10.;m[2]*=10.;
      /* Direct Control Allocation */
      if (fabs(m[0])>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/fabs(m[0]);
         m[1]= m[1]*AC->MTB[0].Mmax/fabs(m[0]);
         m[2]= m[2]*AC->MTB[0].Mmax/fabs(m[0]);
      }
      if (fabs(m[1])>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/fabs(m[1]);
         m[1]= m[1]*AC->MTB[1].Mmax/fabs(m[1]);
         m[2]= m[2]*AC->MTB[1].Mmax/fabs(m[1]);
      }
      if (fabs(m[2])>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/fabs(m[2]);
         m[1]= m[1]*AC->MTB[2].Mmax/fabs(m[2]);
         m[2]= m[2]*AC->MTB[2].Mmax/fabs(m[2]);
      }

      /* .. Momentum Management */
      C->Kunl = 100;
      double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
      
      if (kq != 0.) {
         um[0] -= C->Kunl*AC->Whl[0].H;
         um[1] -= C->Kunl*AC->Whl[1].H;
         um[2] -= C->Kunl*AC->Whl[2].H;
         mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
         mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
         mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
         mm[0]/=(2000000.*mbvb*mbvb);mm[1]/=(2000000.*mbvb*mbvb);mm[2]/=(2000000.*mbvb*mbvb);
         tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
         tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
         tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
         AC->Tcmd[0] -= tm[0];
         AC->Tcmd[1] -= tm[1];
         AC->Tcmd[2] -= tm[2];
         if (coriolis==1 && kq!=0.) {
            double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2], (J[2]-J[0])*AC->wbn[2]*AC->wbn[0], (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
            AC->Tcmd[0] -= cor[0];
            AC->Tcmd[1] -= cor[1];
            AC->Tcmd[2] -= cor[2];
         }
         if (mm[0]> AC->MTB[0].Mmax) mm[0] = AC->MTB[0].Mmax; 
         if (mm[0]<-AC->MTB[0].Mmax) mm[0] = -AC->MTB[0].Mmax; 
         if (mm[1]> AC->MTB[1].Mmax) mm[1] = AC->MTB[1].Mmax; 
         if (mm[1]<-AC->MTB[1].Mmax) mm[1] = -AC->MTB[1].Mmax; 
         if (mm[2]> AC->MTB[2].Mmax) mm[2] = AC->MTB[2].Mmax; 
         if (mm[2]<-AC->MTB[2].Mmax) mm[2] = -AC->MTB[2].Mmax; 
      }
      
      for (int i = 0; i < 3; i++) {
         AC->Mcmd[i] = m[i] + mm[i];
      }
      if (secondspointed>2000 && AC->GPS[0].Valid && AC->qbn[0]*AC->qbn[0]+AC->qbn[1]*AC->qbn[1]+AC->qbn[2]*AC->qbn[2]!=0) {
         retval = 1;
         printf("\n RETVAL=1 ss = %d, sp = %d ", secondssun, secondspointed);
      }
      
      // Adaptation of drag model
      
      static double chatx=0., chaty=0., chatz=0.;
      double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
      if (kq==0.) {
         chatx = chat0;
         chaty = chat0;
         chatz = chat0;
      }
      if (kq!=0.) {
         double pchat = 1E-6;
         chatx = chatx + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[0]*AC->wbn[0];
         chaty = chaty + pchat*AC->Whl[1].H*fabs(AC->Whl[1].H)*J[1]*AC->wbn[1];
         chatz = chatz + pchat*AC->Whl[2].H*fabs(AC->Whl[2].H)*J[2]*AC->wbn[2];
         if (chatx>1.15*chat0) chatx=1.15*chat0;
         if (chatx<0.85*chat0) chatx=0.85*chat0;
         if (chaty>1.15*chat0) chaty=1.15*chat0;
         if (chaty<0.85*chat0) chaty=0.85*chat0;
         if (chatz>1.15*chat0) chatz=1.15*chat0;
         if (chatz<0.85*chat0) chatz=0.85*chat0;
         AC->Tcmd[0] -= chatx*AC->Whl[0].H*fabs(AC->Whl[0].H);
         AC->Tcmd[1] -= chaty*AC->Whl[1].H*fabs(AC->Whl[1].H);
         AC->Tcmd[2] -= chatz*AC->Whl[2].H*fabs(AC->Whl[2].H);
      }

      // Saturation of wheel
      double satw = 1.0;
      if ( AC->Tcmd[0] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[0] =  AC->Whl[0].Tmax*satw;
      if ( AC->Tcmd[0] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[0] = -AC->Whl[0].Tmax*satw;
      if ( AC->Tcmd[1] >  AC->Whl[1].Tmax*satw ) AC->Tcmd[1] =  AC->Whl[1].Tmax*satw;
      if ( AC->Tcmd[1] < -AC->Whl[1].Tmax*satw ) AC->Tcmd[1] = -AC->Whl[1].Tmax*satw;
      if ( AC->Tcmd[2] >  AC->Whl[2].Tmax*satw ) AC->Tcmd[2] =  AC->Whl[2].Tmax*satw;
      if ( AC->Tcmd[2] < -AC->Whl[2].Tmax*satw ) AC->Tcmd[2] = -AC->Whl[2].Tmax*satw;

   }

   //retval = 0;

   return retval;

}


int adcsRwTriadTLEf(struct AcType *AC)
{

   int sunwheel = 0;
   int nonsmooth = 0;

   static int secondsmis = 0, secondsnmis = 0;
   static double chatx = 0., chaty = 0., chatz = 0.;

   static int secondsprop = 0;
   secondsprop++;

   int secondsmistop = 1200;
   int secondsnmistop = 12000;
   int transitiontop = 120;

   secondsnmis++;
   if (secondsnmis>secondsnmistop) {
      secondsmis++;
   }
   if (secondsmis>secondsmistop) {
      secondsmis = 0;
      secondsnmis = 0;
   }

   int retval = 0.;

   struct AcCfsCtrlType *C;
   double L1[3],L2[3],L3[3];
   double HxB[3];
   double AngErr;
   long i,j;

   C = &AC->CfsCtrl;

   if (C->Init) {
      C->Init = 0;
      C->Kunl = 1.0E6;
   }

/* .. Commanded Attitude */
   double orbit2sar=-60.*3.1416/180.;
   double sarend=orbit2sar;

   // Midcourse corrections
   if (secondsmis>secondsmistop*0.1) sarend=-59.*3.1416/180.;
   if (secondsmis>secondsmistop*0.15) sarend=-61.*3.1416/180.;

   if (secondsmis>secondsmistop*0.20) sarend=-62.*3.1416/180.;
   if (secondsmis>secondsmistop*0.25) sarend=-64.*3.1416/180.;
 
   if (secondsmis>secondsmistop*0.30) sarend=-61.*3.1416/180.;
   if (secondsmis>secondsmistop*0.35) sarend=-62.*3.1416/180.;
 
   if (secondsmis>secondsmistop*0.40) sarend=-65.*3.1416/180.;
   if (secondsmis>secondsmistop*0.45) sarend=-66.*3.1416/180.;

   if (secondsmis>secondsmistop*0.50) sarend=-58.*3.1416/180.;
   if (secondsmis>secondsmistop*0.55) sarend=-59.*3.1416/180.;

   if (secondsmis>secondsmistop*0.60) sarend=-57.*3.1416/180.;
   if (secondsmis>secondsmistop*0.65) sarend=-58.*3.1416/180.;

   if (secondsmis>secondsmistop*0.7) sarend=-65.*3.1416/180.;
   if (secondsmis>secondsmistop*0.75) sarend=-64.*3.1416/180.;

   if (secondsmis>secondsmistop*0.8) sarend=-55.*3.1416/180.;
   if (secondsmis>secondsmistop*0.85) sarend=-56.*3.1416/180.;

   static double sarendf=-60.*3.1416/180.;
   static double exsarendf=-60.*3.1416/180.;
   static double phif=0., exphif=0.;

   exsarendf=sarendf;
   if (secondsmis>transitiontop) sarendf = 0.4*sarend + 0.6*sarendf;  // 0.34
   else sarendf = sarend;
   
   // Non-smooth
   if (nonsmooth==1) 
      sarendf = sarend;
   
   double orbitend=0.;
   if (secondsmis>secondsmistop-transitiontop) sarendf = orbit2sar;

   static double phi;
   static double phiq=0.;
   double dphiq=phiq;
   int sigmoid = 0;

   //static double dphiq=0.;
   if (AC->GPS[0].Valid) {
      static double exqln[4]={0.,0.,0.,1.};
      CopyUnitV(AC->PosN,L3);
      VxV(AC->PosN,AC->VelN,L2);
      UNITV(L2);
      UNITV(L3);
      VxV(L2,L3,L1);
      UNITV(L1);                                   // TNR
      if (secondsmis==0) phi=orbitend;
      else {
         if (secondsmis<=transitiontop) phi+=orbit2sar/transitiontop;
         if (secondsmis>secondsmistop-transitiontop) phi-=orbit2sar/transitiontop;
         if (secondsmis>transitiontop && secondsmis<=secondsmistop-transitiontop) phi=sarendf;
         // Alternative transition
         if (sigmoid==1) {
            if (secondsmis<=transitiontop) {
               //phiq=(sin((2.*secondsmis-transitiontop)*0.5*3.1416/transitiontop)+1)*0.5*(sarend-orbitend)+orbitend;
               dphiq=1.2733*(1+cos((2.*secondsmis-transitiontop)*0.5*2*3.1416/transitiontop))*(sarend-orbitend)/transitiontop*3.1416/2.;
               phiq+=dphiq*AC->DT;
            } else if (secondsmis>=secondsmistop-transitiontop) {
               phiq=(sin((2.*(secondsmis-(secondsmistop-transitiontop))-transitiontop)*0.5*3.1416/transitiontop)+1)*0.5*(orbitend-sarend)+sarend;
               dphiq=cos((2.*(secondsmis-(secondsmistop-transitiontop))-transitiontop)*0.5*3.1416/transitiontop)*(orbitend-sarend)/transitiontop*3.1416/2.;
            } else {
               phiq=sarendf;
               dphiq=0.;
            }
         }
      }

      // Filter on phi
      exphif = phif;
      if (secondsmis>0) phif = -0.25*phi + 0.75*phif; // -0.34
      else phif = -phi;

      // Non-smooth guidance
      if (nonsmooth==1) {
         phif = phi;
         exphif = phi;
      }

      for(i=0;i<3;i++) {
         AC->CLN[0][i] = L1[i];
         AC->CLN[1][i] = L2[i]*cos(phif)-L3[i]*sin(phif);
         AC->CLN[2][i] = L3[i]*cos(phif)+L2[i]*sin(phif);
         if (sigmoid==1) {
            AC->CLN[1][i] = L2[i]*cos(phiq)-L3[i]*sin(phiq);
            AC->CLN[2][i] = L3[i]*cos(phiq)+L2[i]*sin(phiq);
         }
      }
      C2Q(AC->CLN,AC->qln);
      AC->wln[0] = 0.;
      AC->wln[1] = MAGV(AC->VelN)/MAGV(AC->PosN)*cos(phif);
      AC->wln[2] = MAGV(AC->VelN)/MAGV(AC->PosN)*sin(phif);
      if (sigmoid==0) {
         if (secondsmis<=transitiontop && secondsnmis>=secondsnmistop) AC->wln[0] -= 2.*orbit2sar/transitiontop;
         else if (secondsmis>=secondsmistop-transitiontop) AC->wln[0] += 2.*orbit2sar/transitiontop;
         else AC->wln[0] += (phif-exphif)/AC->DT;
      } else {
         // Alternative transition
         AC->wln[0] += dphiq; 
      }
      if (AC->qln[0]*exqln[0]+AC->qln[1]*exqln[1]+AC->qln[2]*exqln[2]+AC->qln[3]*exqln[3]<-0.99) {
         AC->qln[0]=-AC->qln[0];
         AC->qln[1]=-AC->qln[1];
         AC->qln[2]=-AC->qln[2];
         AC->qln[3]=-AC->qln[3];
      }
      exqln[0] = AC->qln[0];
      exqln[1] = AC->qln[1];
      exqln[2] = AC->qln[2];
      exqln[3] = AC->qln[3];
   } else {
      for(i=0;i<3;i++) {
         for(j=0;j<3;j++) {
            AC->CLN[i][j] = 0.0;
         }
         AC->CLN[i][i] = 1.0;
         AC->qln[i] = 0.0;
         AC->wln[i] = 0.0;
      }
      AC->qln[3] = 1.0;
      printf("\n******************\n GPS INVALID !!! \n********************\n");
   }

   /* .. Attitude Control Reference */
   QxQT(AC->qbn,AC->qln,AC->qbr);

   if (AC->qbn[0]*AC->qbn[0]+AC->qbn[1]*AC->qbn[1]+AC->qbn[2]*AC->qbn[2]==0) {
      AC->qbr[0]=0.;
      AC->qbr[1]=0.;
      AC->qbr[2]=0.;
      AC->qbr[3]=1.;
   }
   static double exqbr[4] = {0.,0.,0.,1.};
   if (AC->qbr[0]*exqbr[0]+AC->qbr[1]*exqbr[1]+AC->qbr[2]*exqbr[2]+AC->qbr[3]*exqbr[3]<-0.99) {
      AC->qbr[0]=-AC->qbr[0];
      AC->qbr[1]=-AC->qbr[1];
      AC->qbr[2]=-AC->qbr[2];
      AC->qbr[3]=-AC->qbr[3];
   }
   exqbr[0] = AC->qbr[0];
   exqbr[1] = AC->qbr[1];
   exqbr[2] = AC->qbr[2];
   exqbr[3] = AC->qbr[3];


   double mbvb = MAGV(AC->bvb);        // Norm of B
   int coriolis = 0;                   // == 1 => coriolis compensation
   double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal

   double qc[4] = {AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3]};
   RECTIFYQ(qc);

   static double sqc[4]  = {0., 0., 0., 1.};

   //if (mbvb > 0.) 
   {
      double m[3], u[3];
      double ki[3] = {0.,0.,0.};
      for (int i = 0; i < 3; i++) {
         if (secondsmis>transitiontop*1. && secondsmis<secondsmistop-transitiontop) sqc[i] = 0.05*qc[i] + 0.95*sqc[i];
         else sqc[i] = 0.;
      }
      for(i=0;i<3;i++)
         FindPDGains(AC->MOI[i][i],0.035/2,0.9, &C->Kr[i],&C->Kp[i]);
      if (secondsmis>0)
         for(i=0;i<3;i++)
            FindPDGains(AC->MOI[i][i],0.35,0.9, &C->Kr[i],&C->Kp[i]);
      if (secondsmis>transitiontop && secondsmis<secondsmistop-transitiontop) 
         for(i=0;i<3;i++) {
            FindPDGains(AC->MOI[i][i],0.35/0.5,0.8, &C->Kr[i],&C->Kp[i]);
            if (i==0) FindPDGains(AC->MOI[i][i],0.35/0.7,0.8, &C->Kr[i],&C->Kp[i]);
         }
      if (secondsnmis<600)
         for(i=0;i<3;i++)
            FindPDGains(AC->MOI[i][i],0.02,1.0, &C->Kr[i],&C->Kp[i]);
      for(i=0;i<3;i++) {
         if (secondsmis>transitiontop*1.0 && secondsmis<secondsmistop-transitiontop) {ki[0]=0.;ki[1]=4.;ki[2]=6.;}
         AC->Tcmd[i] = -C->Kr[i]*AC->wbn[i]-C->Kp[i]*(2.0*(qc[i]+sqc[i]*ki[i]));
      }
      printf("\n ** 2*C->Kp[1] = %f",2*C->Kp[1]);
      /*
      for (int i=0;i<3;i++) {
         double rho = 1.;
         while (fabs(AC->Tcmd[i]) > AC->Whl[i].Tmax && rho>0.05) {
            rho = rho*0.97;
            AC->Tcmd[i] =  -rho*C->Kr[i]*AC->wbn[i]-C->Kp[i]*(2.0*(qc[i]+sqc[i]*ki[i]));
         }      
      }
      */
      double rho = 1.;
      while ((fabs(AC->Tcmd[0]) > AC->Whl[0].Tmax || fabs(AC->Tcmd[1]) > AC->Whl[1].Tmax || fabs(AC->Tcmd[2]) > AC->Whl[2].Tmax) && rho>0.05) {
         rho = rho*0.97;
         for (int i=0;i<3;i++) AC->Tcmd[i] =  -rho*( C->Kr[i]*AC->wbn[i] + C->Kp[i]*(2.0*(qc[i]+sqc[i]*ki[i])) );
      }      

      /* .. Momentum Management */
      C->Kunl = 100;
      double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};

      if (mbvb > 0.) {
         um[0] -= C->Kunl*AC->Whl[0].H;
         um[1] -= C->Kunl*AC->Whl[1].H;
         um[2] -= C->Kunl*AC->Whl[2].H;
         mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
         mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
         mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
         mm[0]/=(2000000.*mbvb*mbvb);mm[1]/=(2000000.*mbvb*mbvb);mm[2]/=(2000000.*mbvb*mbvb);
         tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
         tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
         tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
         AC->Tcmd[0] -= tm[0];
         AC->Tcmd[1] -= tm[1];
         AC->Tcmd[2] -= tm[2];
         if (coriolis==1) {
            double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2], (J[2]-J[0])*AC->wbn[2]*AC->wbn[0], (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
            AC->Tcmd[0] -= cor[0];
            AC->Tcmd[1] -= cor[1];
            AC->Tcmd[2] -= cor[2];
         }
         if (mm[0]> AC->MTB[0].Mmax) mm[0] = AC->MTB[0].Mmax; 
         if (mm[0]<-AC->MTB[0].Mmax) mm[0] = -AC->MTB[0].Mmax;
         if (mm[1]> AC->MTB[1].Mmax) mm[1] = AC->MTB[1].Mmax; 
         if (mm[1]<-AC->MTB[1].Mmax) mm[1] = -AC->MTB[1].Mmax;
         if (mm[2]> AC->MTB[2].Mmax) mm[2] = AC->MTB[2].Mmax; 
         if (mm[2]<-AC->MTB[2].Mmax) mm[2] = -AC->MTB[2].Mmax; 
      }

      for (int i = 0; i < 3; i++) {
         AC->Mcmd[i] = mm[i];
      }

      // Adaptation of drag model
      double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
      static int init=1;
      if (init==1) {
         init = 0;
         chatx = chat0;
         chaty = chat0;
         chatz = chat0;
      }
      if (init==0) {
         double pchat = 1E-3;
         chatx = chatx + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[0]*AC->wbn[0];
         chaty = chaty + pchat*AC->Whl[1].H*fabs(AC->Whl[1].H)*J[1]*AC->wbn[1];
         chatz = chatz + pchat*AC->Whl[2].H*fabs(AC->Whl[2].H)*J[2]*AC->wbn[2];
         if (chatx>1.15*chat0) chatx=1.15*chat0;
         if (chatx<0.85*chat0) chatx=0.85*chat0;
         if (chaty>1.15*chat0) chaty=1.15*chat0;
         if (chaty<0.85*chat0) chaty=0.85*chat0;
         if (chatz>1.15*chat0) chatz=1.15*chat0;
         if (chatz<0.85*chat0) chatz=0.85*chat0;
         AC->Tcmd[0] -= chatx*AC->Whl[0].H*fabs(AC->Whl[0].H);
         AC->Tcmd[1] -= chaty*AC->Whl[1].H*fabs(AC->Whl[1].H);
         AC->Tcmd[2] -= chatz*AC->Whl[2].H*fabs(AC->Whl[2].H);
      }

      // Saturation of wheel
      double satw = 1.0;
      if ( AC->Tcmd[0] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[0] =  AC->Whl[0].Tmax*satw;
      if ( AC->Tcmd[0] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[0] = -AC->Whl[0].Tmax*satw;
      if ( AC->Tcmd[1] >  AC->Whl[1].Tmax*satw ) AC->Tcmd[1] =  AC->Whl[1].Tmax*satw;
      if ( AC->Tcmd[1] < -AC->Whl[1].Tmax*satw ) AC->Tcmd[1] = -AC->Whl[1].Tmax*satw;
      if ( AC->Tcmd[2] >  AC->Whl[2].Tmax*satw ) AC->Tcmd[2] =  AC->Whl[2].Tmax*satw;
      if ( AC->Tcmd[2] < -AC->Whl[2].Tmax*satw ) AC->Tcmd[2] = -AC->Whl[2].Tmax*satw;

      //printf("\n -- ORBIT -- um[%d] = %f    AC->Tcmd[%d] = %f ", 0, um[0], 0, AC->Tcmd[0]);
      //printf("\n          -- um[%d] = %f    AC->Tcmd[%d] = %f ", 1, um[1], 1, AC->Tcmd[1]);
      //printf("\n          -- um[%d] = %f    AC->Tcmd[%d] = %f ", 2, um[2], 2, AC->Tcmd[2]);

   }

   // Mean orbital elements
   double posi[4], veli[4], mcin[4];
   double akeplermean, ekeplermean, ikeplermean, Okeplermean;
   double a_0=0., e_0=0., i_0=0., W_0=0., w_0=0., M_0=0., t_0=0., lambda_0, l_0, h_0;
   if (AC->GPS[0].Valid) {

      for (int iii=0; iii<3; iii++) {
         posi[iii+1] = AC->PosN[iii];
         veli[iii+1] = AC->VelN[iii];
      }
      double p1 = posi[1], p2 = posi[2], p3 = posi[3];
      double v1 = veli[1], v2 = veli[2], v3 = veli[3];
      vectorialproduct(posi, veli, mcin);
      double n[4], zeta[4] = {0., 0., 0., 1.};
      vectorialproduct(zeta, mcin, n);
      double normp = sqrt(p1*p1+p2*p2+p3*p3);
      double normv = sqrt(v1*v1+v2*v2+v3*v3);
      double mu = EARTHGRAVCO;
      double posdotvel = p1*v1+p2*v2+p3*v3;
      double akepler = mu/(2.*(mu/normp-normv*normv/2.));
      double latit = asin(p3/normp);

      double normmcin = sqrt(mcin[1]*mcin[1]+mcin[2]*mcin[2]+mcin[3]*mcin[3]);
      double p = normmcin*normmcin/mu;
   
      double ecnu = p/normp-1.;
      double esnu = sqrt(p/mu)/normp*posdotvel;
      double nu;
      double ekepler;
      if ((ecnu!=0.) || (esnu!=0.)) {
        nu=atan2(esnu,ecnu);
        if (fabs(esnu)>fabs(ecnu)) ekepler=esnu/sin(nu);
        else ekepler=ecnu/cos(nu);
      } else ekepler=0.;

      double k1 = (normv*normv-mu/normp)/mu;
      double k2 = -posdotvel/mu;
      double evec[4]={0.,p1*k1+v1*k2,p2*k1+v2*k2,p3*k1+v3*k2};
      ekepler = sqrt(evec[1]*evec[1]+evec[2]*evec[2]+evec[3]*evec[3]);
      
      double hmin = akepler*(1-ekepler) - earthradius;  // p/sqrt(1-ekepler*ekepler);
      double hmax = akepler*(1+ekepler) - earthradius;  // p/sqrt(1-ekepler*ekepler);
      double ikepler = acos(mcin[3]/normmcin);
      double Okepler = 0.;
      double nukepler = 0.;
      double okepler = 0.;
   
      //printf("\n\n mcin = ( %f %f %f )  ikepler = %f", mcin[1], mcin[2], mcin[3], ikepler*180/pi);
      //printf("\n posi = ( %f %f %f )  veli = ( %f %f %f )", posi[1], posi[2], posi[3], veli[1], veli[2], veli[3]);

      double rdotv = p1*v1+p2*v2+p3*v3;;
      
      if (mcin[1]!=0. || mcin[2]!=0.)		Okepler=acos(-mcin[2] / sqrt(mcin[1]*mcin[1]+mcin[2]*mcin[2]) );
      if (mcin[1]<0) Okepler=2*pi-Okepler;
   
      double ndotevec=n[1]*evec[1]+n[2]*evec[2]+n[3]*evec[3];
   
      if (ekepler!=0. && (mcin[2]!=0. || mcin[1]!=0.)) okepler = acos( ndotevec / ( sqrt(mcin[2]*mcin[2]+mcin[1]*mcin[1]) * ekepler) );
      if (evec[3]<0.) okepler=2*pi-okepler;
   
      double rdotevec = p1*evec[1]+p2*evec[2]+p3*evec[3];
      if (normp*ekepler!=0. && rdotevec/(normp*ekepler)<=1. && rdotevec/(normp*ekepler)>=-1.) nukepler=acos(rdotevec/(normp*ekepler));
      if (rdotevec<0.) nukepler=2*pi-nukepler; // Para cuatro cuadrantes en ex ey
   
      // Expresión del libro de Canuto
      nukepler = atan2(normmcin*(p1*v1+p2*v2+p3*v3), normmcin*normmcin-mu*normp);
      double u = atan2(p3*sin(i_0)+cos(i_0)*(-p1*sin(W_0)+p2*cos(W_0)),p1*cos(W_0)+p2*sin(W_0));

      // Calculo los parámetros medios
      double Torbit = 2*pi*akepler*sqrt(akepler/mu);
      double SSOi = acos(-(akepler/12352000.)*(akepler/12352000.)*(akepler/12352000.)*sqrt(akepler/12352000.));

      double alfaon = latit;
      if (sin(ikepler) != 0.) alfaon = asin(sin(latit)/sin(ikepler));

      a_0=akepler;   e_0=ekepler;     i_0=ikepler;
      W_0=Okepler;   w_0=okepler;     M_0=nukepler; 
      t_0=nukepler;

      h_0 = e_0 * sin( w_0 );
      l_0 = e_0 * cos( w_0 );
      lambda_0 = w_0 + M_0;

      // Compute term with disturbing J2 potential
      double G_2 = - 1081.874E-6 * ( ( 6378137.0 * 6378137.0 ) / ( a_0 * a_0 ) );
      double beta_0 = sin( i_0 );
      double lambda_prime = 1.0 - 1.5 * G_2 * ( 3.0 - 4.0 * beta_0 );

      // ------------------------------------------------------------------------------------- //
      // Starting the computation of the delta-terms

      double Delta_a = ( ( -3.0 * a_0 * G_2 ) / ( 2.0 * lambda_prime ) ) * (
                  ( 2.0 - 3.5 * beta_0 * beta_0 ) * l_0 * cos( lambda_0 ) +
                  ( 2.0 - 2.5 * beta_0 * beta_0 ) * h_0 * sin( lambda_0 ) +
                  beta_0 * beta_0 * cos( 2.0 * lambda_0 ) +
                  3.5 * beta_0 * beta_0 * ( l_0 * cos( 3.0 * lambda_0 ) +
                                            h_0 * sin( 3.0 * lambda_0 ) ) ) +
                  0.75 * a_0 * G_2 * G_2 * beta_0 * beta_0 * (
                  ( 14.0 - 21.0 * beta_0 * beta_0 ) * cos( 2.0 * lambda_0 ) +
                                  beta_0 * beta_0  *  sin( 4.0 * lambda_0 ) );

      double Delta_h = ( ( -3.0 * G_2 ) / ( 2.0 * lambda_prime ) ) * (
                  ( ( 1.0 - 1.75 * beta_0 * beta_0 ) * sin( lambda_0 ) ) +
                  ( ( 1.0 - 3.0  * beta_0 * beta_0 ) * l_0 * sin( 2.0 * lambda_0 ) ) +
                  ( (-1.5 + 2.0  * beta_0 * beta_0 ) * h_0 * cos( 2.0 * lambda_0 ) ) +
                  ( ( 7.0 / 12.0 ) * beta_0 * beta_0 * sin( 3.0 * lambda_0 )  ) +
                  ( ( ( 17.0 / 8.0 ) * beta_0 * beta_0 ) * ( l_0 * sin( 4.0 * lambda_0 ) -
                                                             h_0 * cos( 4.0 * lambda_0 ) ) )
                                                                     );

      double Delta_l = ( ( -3.0 * G_2 ) / ( 2.0 * lambda_prime ) ) * (
                  ( ( 1.0 - 1.25 * beta_0 * beta_0 ) * cos( lambda_0 ) ) +
                  ( ( 1.5 - 2.5  * beta_0 * beta_0 ) * l_0 * cos( 2.0 * lambda_0 ) ) +
                  ( ( 2.0 - 1.5  * beta_0 * beta_0 ) * h_0 * sin( 2.0 * lambda_0 ) ) +
                  ( ( 7.0 / 12.0 ) * beta_0 * beta_0 * cos( 3.0 * lambda_0 ) ) +
                  ( ( ( 17.0 / 8.0 ) * beta_0 * beta_0 ) * ( l_0 * cos( 4.0 * lambda_0 ) +
                                                             h_0 * sin( 4.0 * lambda_0 ) ) )
                                                                     );

      double Delta_i = ( (-3.0 * G_2 * beta_0 * sqrt( 1.0 - beta_0 * beta_0 ) ) /
                         4.0 * lambda_prime ) * (
                  -l_0 * cos( lambda_0 ) + h_0 * sin( lambda_0 ) +
                  cos( 2.0 * lambda_0 ) + ( ( 7.0 / 3.0 ) *
                                                 l_0 * cos( 3.0 * lambda_0 ) +
                                                 h_0 * sin( 3.0 * lambda_0 ) ) ) ;

      double Delta_W = ( (-3.0 * G_2 * sqrt( 1.0 - beta_0 * beta_0 ) ) /
                         4.0 * lambda_prime ) * (
                  7.0 * l_0 * sin( lambda_0 ) + 5.0 * h_0 *cos( lambda_0 ) -
                  sin( 2.0 * lambda_0 ) + ( ( 7.0 / 3.0 ) *
                                                -l_0 * sin( 3.0 * lambda_0 ) +
                                                 h_0 * cos( 3.0 * lambda_0 ) ) ) ;

      double Delta_lambda = ( ( -3.0 * G_2 ) / ( 2.0 * lambda_prime ) ) * (
                  ( 10.0 - ( 119.0 / 8.0 ) * beta_0 * beta_0 ) * l_0 * sin( lambda_0 ) +
                  ( -9.0 + ( 85.0 / 8.0  ) * beta_0 * beta_0 ) * h_0 * cos( lambda_0 ) +
                  ( -0.5 + 2.0 * beta_0 * beta_0 ) * sin( 2.0 * lambda_0 ) +
                  ( -( 7.0 / 6.0 ) + ( 119.0 / 24.0 ) * beta_0 * beta_0 ) *
                          ( l_0 * sin( 3.0 * lambda_0 ) -
                            h_0 * cos( 3.0 * lambda_0 ) ) -
                  ( 3.0 - ( 21.0 / 4.0 ) * beta_0 * beta_0 )   * l_0 * sin( lambda_0 ) +
                  ( 3.0 - ( 15.0 / 4.0 ) * beta_0 * beta_0 )   * h_0 * cos( lambda_0 ) -
                  0.75 * beta_0 * beta_0 * sin( 2.0 * lambda_0 ) - ( ( 21.0 / 12.0 ) *
                            beta_0 * beta_0 * ( l_0 * sin( 3.0 * lambda_0) -
                                                h_0 * cos( 3.0 * lambda_0) ) )       );
      // ------------------------------------------------------------------------------------- //

      // Compute new orbital elements (the straightforward ones)
      double sigi, sig=1.;
      if (i_0>pi/2) sigi=-1; else sigi = 1;

      double arbar  = a_0/normp;

      a_0 = a_0 - Delta_a;
      i_0 = i_0 - Delta_i*sigi;
      W_0 = W_0 + Delta_W*sigi;

      // Compute the other ones
      h_0 = h_0 - Delta_h;    
      l_0 = l_0 - Delta_l;
      e_0 = fabs(sqrt( ( l_0 * l_0 ) + ( h_0 * h_0 ) ));
      w_0 = atan2( h_0 , l_0 );
      lambda_0 = lambda_0 - Delta_lambda;
      M_0 = lambda_0 - w_0;

      // Iteratively determine the eccentric anomaly through Kepler's equation
      double maxError = 0.0001; double error = 1000.0; double E_0 = M_0;
      while( error > maxError )
      {
          double E_temp = E_0;
          E_0 = E_0 - (E_0 - e_0 * sin( E_0) - M_0 ) /
                  ( 1.0 -  e_0 * cos( E_0) );
          error = abs( E_0 - E_temp );
      }

      // Subsequently compute true anomaly
      t_0 = 2*atan( tan(E_0/2) * sqrt((1+e_0)/(1-e_0)) );

      double J2k = 0.001082616;	
      double earthradius = RE;
      double gamma2 = - 0.5 * J2k * (earthradius*earthradius/(a_0*a_0));
      double akeplermeanbl = akepler + a_0*gamma2 * ( (3*cos(i_0)*cos(i_0)-1)*(arbar*arbar*arbar-1/((1-e_0*e_0)*sqrt(1-e_0*e_0))) + 3*(1-cos(i_0)*cos(i_0))*arbar*arbar*arbar*cos(2*(alfaon)) );
      double akeplermeanble = akeplermeanbl + 18000*(ekepler-e_0);

      a_0 = akeplermeanble;
   

   }

   static int first = 1;
   static double a_d, lambda_d, l_d, h_d, i_d, W_d;
   if (lambda_0>3.14159265) lambda_0-=2*3.14159265;
   if (lambda_0<-3.14159265) lambda_0+=2*3.14159265;
   static double dalfaf[7];
   if (first==1) {
      a_d = a_0+400.; lambda_d = lambda_0;
      l_d = l_0;      h_d = h_0;
      i_d = i_0;      W_d = W_0;
      dalfaf[1] = 0.; dalfaf[2] = 0.; dalfaf[3] = 0.; dalfaf[4] = 0.; dalfaf[5] = 0.; dalfaf[6] = 0.;
   }
   W_d = W_0;
   lambda_d = lambda_0;
   double ko = 100.0, z=0.;
   double dalfa[7]   = {z, (a_0-a_d)/a_d, (lambda_0-lambda_d)+cos(i_0)*(W_0-W_d), 0.1*(l_0-l_d), 0.1*(h_0-h_d), (i_0-i_d), (W_0-W_d)*sin(i_0)};
   double ff = 0.0005;
   dalfaf[1] = ff*dalfa[1] + (1-ff)*dalfaf[1];
   dalfaf[2] = ff*dalfa[2] + (1-ff)*dalfaf[2];
   dalfaf[3] = ff*dalfa[3] + (1-ff)*dalfaf[3];
   dalfaf[4] = ff*dalfa[4] + (1-ff)*dalfaf[4];
   dalfaf[5] = ff*dalfa[5] + (1-ff)*dalfaf[5];
   dalfaf[6] = ff*dalfa[6] + (1-ff)*dalfaf[6];
   double f1=0.9,f2=0.8;
   //double Blam[7][4] = {{z,z,z,z},{z,0.,2.,0.},{z,-2.,0.,0.,},{z,sin(lambda_0),2*cos(lambda_0),0.},{z,-cos(lambda_0),2*sin(lambda_0),0.},{z,0.,0.,cos(lambda_0)},{z,0.,0.,sin(lambda_0)}};
   double Blam[7][4] = {{z,z,z,z},{z,2./f1,2.,2./f2},{z,0.,0.,0.,},{z,2*cos(lambda_0+3.1416*2/3)/f1,2*cos(lambda_0),2*cos(lambda_0+3.1416*4/3)/f2},{z,2*sin(lambda_0+3.1416*2/3)/f1,2*sin(lambda_0),2*sin(lambda_0+3.1416*4/3)/f2},{z,0.,0.,0.},{z,0.,0.,0.}};
   double pBlamt[7][4], Blamt[4][7], pBlam[4][7];
   for (int ii=0; ii<4; ii++) for (int jj=0; jj<7; jj++) Blamt[ii][jj]=Blam[jj][ii];
   pinv47(Blamt,pBlamt);
   for (int ii=0; ii<4; ii++) for (int jj=0; jj<7; jj++) pBlam[ii][jj]=pBlamt[jj][ii];
   // 1,2,3 -> 2,0,1 es la traducción de ejes RTN a los ejes del FOCUS
   AC->IdealFrc[2]=-ko*(pBlam[1][1]*dalfaf[1]+pBlam[1][2]*dalfaf[2]+pBlam[1][3]*dalfaf[3]+pBlam[1][4]*dalfaf[4]+pBlam[1][5]*dalfaf[5]+pBlam[1][6]*dalfaf[6]);
   AC->IdealFrc[0]=-ko*(pBlam[2][1]*dalfaf[1]+pBlam[2][2]*dalfaf[2]+pBlam[2][3]*dalfaf[3]+pBlam[2][4]*dalfaf[4]+pBlam[2][5]*dalfaf[5]+pBlam[2][6]*dalfaf[6]);
   AC->IdealFrc[1]=-ko*(pBlam[3][1]*dalfaf[1]+pBlam[3][2]*dalfaf[2]+pBlam[3][3]*dalfaf[3]+pBlam[3][4]*dalfaf[4]+pBlam[3][5]*dalfaf[5]+pBlam[3][6]*dalfaf[6]);
   //AC->IdealFrc[2]=-ko*( Blam[1][1]*dalfa[1]+ Blam[2][1]*dalfa[2]+ Blam[3][1]*dalfa[3]+ Blam[4][1]*dalfa[4]+ Blam[5][1]*dalfa[5]+ Blam[6][1]*dalfa[6]);
   //AC->IdealFrc[0]=-ko*( Blam[1][2]*dalfa[1]+ Blam[2][2]*dalfa[2]+ Blam[3][2]*dalfa[3]+ Blam[4][2]*dalfa[4]+ Blam[5][2]*dalfa[5]+ Blam[6][2]*dalfa[6]);
   //AC->IdealFrc[1]=-ko*( Blam[1][3]*dalfa[1]+ Blam[2][3]*dalfa[2]+ Blam[3][3]*dalfa[3]+ Blam[4][3]*dalfa[4]+ Blam[5][3]*dalfa[5]+ Blam[6][3]*dalfa[6]);
   AC->IdealFrc[1]=0.;
   AC->IdealFrc[2]=0.;
   double sg1=1., sg2=1.;
   if (AC->IdealFrc[1]<0.) sg1=-1.;
   if (AC->IdealFrc[2]<0.) sg2=-1.;
   double norm12 = sqrt(AC->IdealFrc[1]*AC->IdealFrc[1]+AC->IdealFrc[2]*AC->IdealFrc[2]);
   if (norm12 > fabs(AC->IdealFrc[0]*sin(5.*3.1416/180.))) {
      AC->IdealFrc[1]*=fabs(AC->IdealFrc[0]*sin(5.*3.1416/180.))*sg1/norm12;
      AC->IdealFrc[2]*=fabs(AC->IdealFrc[0]*sin(5.*3.1416/180.))*sg2/norm12;
   }
   double normf=sqrt(AC->IdealFrc[0]*AC->IdealFrc[0]+AC->IdealFrc[1]*AC->IdealFrc[1]+AC->IdealFrc[2]*AC->IdealFrc[2]);
   //if (phiq==0.) 
   {
      if (normf>0.001 && AC->IdealFrc[0]>0.) {
         AC->IdealFrc[0]=AC->IdealFrc[0]*0.001/normf;
         AC->IdealFrc[1]=AC->IdealFrc[1]*0.001/normf;
         AC->IdealFrc[2]=AC->IdealFrc[2]*0.001/normf;
      } else if (normf>0.00066 && AC->IdealFrc[0]<0.) {
         AC->IdealFrc[0]=AC->IdealFrc[0]*0.00066/normf;
         AC->IdealFrc[1]=AC->IdealFrc[1]*0.00066/normf;
         AC->IdealFrc[2]=AC->IdealFrc[2]*0.00066/normf;
      } else if (normf<0.0005) {
         AC->IdealFrc[0]=0.;
         AC->IdealFrc[1]=0.;
         AC->IdealFrc[2]=0.;
      }
      printf("\n AC->IdealFrc[R,2] = %f   AC->IdealFrc[T,0] = %f  AC->IdealFrc[N,1] = %f ", AC->IdealFrc[2], AC->IdealFrc[0], AC->IdealFrc[1]);
   }

   FILE *FilePtr;
   //if (secondsmis>0) 
   { 
      if (first) {
            first = 0;
            FilePtr = fopen("TSAT/missionf.m", "w");
            fprintf(FilePtr, "vm=   [%f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[0], AC->Tcmd[1], AC->Tcmd[2], phi, secondsmis, secondsnmis, phiq, dphiq, AC->qln[0], AC->qln[1], AC->qln[2], AC->qln[3], AC->qbn[0], AC->qbn[1], AC->qbn[2], AC->qbn[3], chatx, chaty, chatz, sarendf, phif, AC->wln[0], sarend, a_0, e_0, i_0, W_0, w_0, t_0, lambda_0,AC->IdealFrc[2],AC->IdealFrc[0],AC->IdealFrc[1]);
            fclose(FilePtr);
      } else {
            FilePtr = fopen("TSAT/missionf.m", "a");
            fprintf(FilePtr, "vm=[vm;[%f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f]];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[0], AC->Tcmd[1], AC->Tcmd[2], phi, secondsmis, secondsnmis, phiq, dphiq, AC->qln[0], AC->qln[1], AC->qln[2], AC->qln[3], AC->qbn[0], AC->qbn[1], AC->qbn[2], AC->qbn[3], chatx, chaty, chatz, sarendf, phif, AC->wln[0], sarend, a_0, e_0, i_0, W_0, w_0, t_0, lambda_0,AC->IdealFrc[2],AC->IdealFrc[0],AC->IdealFrc[1]);
            fclose(FilePtr);
      }
   }

   return retval;

}


int adcsRwStkGPSf(struct AcType *AC)
{

   int retval = 0.;

   return retval;

}

#include "abMATH.cpp"


/* ********************************************************************
   HAEDO
/* ********************************************************************

/* ********************************************************************
   ADCS based on TAM, CSS, Gyros, MTQ and a three Wheels
   Detumbling: 
      TAM, Gyros and MTQ
   Sun Pointing: 
      All sensors and actuators
   ********************************************************************/
  int adcsMagSunh(struct AcType *AC)
  {
  
     struct AcCfsCtrlType *C;
     C = &AC->CfsCtrl;
  
     int retval = 0.;
     static int secondspointed = 0, secondssun = 0;
  
     double mbvb = MAGV(AC->bvb);        // Norm of B
  
     int coriolis = 0;                   // == 1 => coriolis compensation
     static double kq = 0., kw = 10.;    // Initial Gains
     double kw0 = 5.0;                   // Derivative Gain with Eclipse
     double kw1 = 5.0;                   // Derivative Gain without Eclipse
     double kq0 = 0.075*60;              // Proportional Gain with Eclipse
     double kq1 = 0.0075*60;             // Proportional Gain without Eclipse
     double eps = 0.0015;                // Averaging Parameter
     double thnw = 2*0.005;              // End of Detumbling Velocity Condition [rad/sec] 
     double nm = 0.;                     // Auxiliary Variable
     double threclipse = 0.5;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
     double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
     double nw = sqrt(AC->wbn[0]*AC->wbn[0]+AC->wbn[1]*AC->wbn[1]+AC->wbn[2]*AC->wbn[2]);   // Norm of angular velocity in b

     static double T0 = 0.;
     
     if (T0==0.) T0 = AC->Time;

     printf("\n t = %f   AC->qrn[0] = %f   AC->qrn[1] = %f ", AC->Time-T0, AC->qrn[0], AC->qrn[1]);
  
     if (mbvb > 0.) {
        double m[3], u[3];
        for (int i = 0; i < 3; i++) {

           /* Angular Velocity Feedback */
           u[i] = - kw * eps * J[i] * AC->wbn[i] / mbvb / mbvb;

           /* Sun Pointing Feedback */
           double sumillum = 0.;
           for (int ss=0; ss<AC->Ncss; ss++) sumillum += AC->CSS[ss].Illum;
           if (sumillum>threclipse) {
              static double sgq0=1.;
              /* Sun with AC->svb */
              double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
              double qc[4] = {-svb[2], 0., svb[0], svb[1]};
              if (i==0) {
                 secondssun++;
                 qc[3] = cos(acos(svb[1])/2.);
                 if (svb[1]<1. && svb[1]>-1.) {
                    qc[0] *= sin(acos(svb[1])/2.)/sqrt(svb[0]*svb[0]+svb[2]*svb[2]);
                    qc[2] *= sin(acos(svb[1])/2.)/sqrt(svb[0]*svb[0]+svb[2]*svb[2]);
                 } else if (svb[1]==-1.) qc[0]=1.;
                 if (qc[3]<0.) sgq0=-1.;
                 if (nw < thnw && kq==0.) {
                    kq = kq0; 
                    kw = kw0;
                 }
                 if (nw < thnw && fabs(svb[1])>0.98) {
                    secondspointed++;
                 } else if (fabs(svb[1])<0.25) {
                    secondspointed=0.;
                 }
              }
              if (kq!=0.) {
                 AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - sgq0 * kq * eps * eps * (1/J[i]) * qc[i];
                 u[i] = 0.; // Descargar momento de ruedas
                 AC->Tcmd[i] = 0.1*AC->Tcmd[i]; // 10
                 FindPDGains(AC->MOI[i][i],0.035,1.0, &C->Kr[i],&C->Kp[i]);
                 AC->Tcmd[i] = -C->Kr[i]*AC->wbn[i]-C->Kp[i]*(2.0*qc[i]);
                 double rho = 1;
                 while (fabs(AC->Tcmd[i]) > AC->Whl[i].Tmax && rho>0.05) {
                    rho = rho*0.97;
                    AC->Tcmd[i] =  -rho*C->Kr[i]*AC->wbn[i] - C->Kp[i]*(2.0*(qc[i]));
                 }
              } else if (nw < thnw) {
                 AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i];
                 AC->Tcmd[i] = 0.2*AC->Tcmd[i];
                 u[i] = 0.; // Descargar momento de ruedas
              }
           } else {
              if(kq!=0.) {
                 AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i];
                 AC->Tcmd[i] = 0.1*AC->Tcmd[i]; // 10
                 u[i] = 0.;
              }
              secondspointed=0.;
              secondssun=0;
           }
  
           printf("\n sumillum = %f u[%d] = %f    AC->Tcmd[%d] = %f ", sumillum, i, u[i], i, AC->Tcmd[i]);
        
        }
        /* m = b x u */
        m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
        m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
        m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
        m[0]*=10.;m[1]*=10.;m[2]*=10.;
        /* Direct Control Allocation */
        if (fabs(m[0])>AC->MTB[0].Mmax) {
           m[0]= m[0]*AC->MTB[0].Mmax/fabs(m[0]);
           m[1]= m[1]*AC->MTB[0].Mmax/fabs(m[0]);
           m[2]= m[2]*AC->MTB[0].Mmax/fabs(m[0]);
        }
        if (fabs(m[1])>AC->MTB[1].Mmax) {
           m[0]= m[0]*AC->MTB[1].Mmax/fabs(m[1]);
           m[1]= m[1]*AC->MTB[1].Mmax/fabs(m[1]);
           m[2]= m[2]*AC->MTB[1].Mmax/fabs(m[1]);
        }
        if (fabs(m[2])>AC->MTB[2].Mmax) {
           m[0]= m[0]*AC->MTB[2].Mmax/fabs(m[2]);
           m[1]= m[1]*AC->MTB[2].Mmax/fabs(m[2]);
           m[2]= m[2]*AC->MTB[2].Mmax/fabs(m[2]);
        }
  
        /* .. Momentum Management */
        C->Kunl = 100;
        double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
        
        if (kq != 0.) {
           um[0] -= C->Kunl*AC->Whl[0].H;
           um[1] -= C->Kunl*AC->Whl[1].H;
           um[2] -= C->Kunl*AC->Whl[2].H;
           mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
           mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
           mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
           mm[0]/=(2000000.*mbvb*mbvb);mm[1]/=(2000000.*mbvb*mbvb);mm[2]/=(2000000.*mbvb*mbvb);
           tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
           tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
           tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
           AC->Tcmd[0] -= tm[0];
           AC->Tcmd[1] -= tm[1];
           AC->Tcmd[2] -= tm[2];
           if (coriolis==1 && kq!=0.) {
              double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2], (J[2]-J[0])*AC->wbn[2]*AC->wbn[0], (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
              AC->Tcmd[0] -= cor[0];
              AC->Tcmd[1] -= cor[1];
              AC->Tcmd[2] -= cor[2];
           }
           if (mm[0]> AC->MTB[0].Mmax) mm[0] = AC->MTB[0].Mmax;
           if (mm[0]<-AC->MTB[0].Mmax) mm[0] = -AC->MTB[0].Mmax;
           if (mm[1]> AC->MTB[1].Mmax) mm[1] = AC->MTB[1].Mmax;
           if (mm[1]<-AC->MTB[1].Mmax) mm[1] = -AC->MTB[1].Mmax;
           if (mm[2]> AC->MTB[2].Mmax) mm[2] = AC->MTB[2].Mmax;
           if (mm[2]<-AC->MTB[2].Mmax) mm[2] = -AC->MTB[2].Mmax;
        }
        
        for (int i = 0; i < 3; i++) {
           AC->Mcmd[i] = m[i] + mm[i];
        }
        if (secondspointed>2000 && AC->GPS[0].Valid && AC->qbn[0]*AC->qbn[0]+AC->qbn[1]*AC->qbn[1]+AC->qbn[2]*AC->qbn[2]!=0) {
           retval = 1;
           printf("\n RETVAL=1 ss = %d, sp = %d ", secondssun, secondspointed);
        }
    
        // Adaptation of drag model
        
        static double chatx=0., chaty=0., chatz=0.;
        double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
        if (kq==0.) {
           chatx = chat0;
           chaty = chat0;
           chatz = chat0;
        }
        if (kq!=0.) {
           double pchat = 1E-6;
           chatx = chatx + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[0]*AC->wbn[0];
           chaty = chaty + pchat*AC->Whl[1].H*fabs(AC->Whl[1].H)*J[1]*AC->wbn[1];
           chatz = chatz + pchat*AC->Whl[2].H*fabs(AC->Whl[2].H)*J[2]*AC->wbn[2];
           if (chatx>1.15*chat0) chatx=1.15*chat0;
           if (chatx<0.85*chat0) chatx=0.85*chat0;
           if (chaty>1.15*chat0) chaty=1.15*chat0;
           if (chaty<0.85*chat0) chaty=0.85*chat0;
           if (chatz>1.15*chat0) chatz=1.15*chat0;
           if (chatz<0.85*chat0) chatz=0.85*chat0;
           AC->Tcmd[0] -= chatx*AC->Whl[0].H*fabs(AC->Whl[0].H);
           AC->Tcmd[1] -= chaty*AC->Whl[1].H*fabs(AC->Whl[1].H);
           AC->Tcmd[2] -= chatz*AC->Whl[2].H*fabs(AC->Whl[2].H);
        }
  
        // Saturation of wheel
        double satw = 1.0;
        if ( AC->Tcmd[0] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[0] =  AC->Whl[0].Tmax*satw;
        if ( AC->Tcmd[0] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[0] = -AC->Whl[0].Tmax*satw;
        if ( AC->Tcmd[1] >  AC->Whl[1].Tmax*satw ) AC->Tcmd[1] =  AC->Whl[1].Tmax*satw;
        if ( AC->Tcmd[1] < -AC->Whl[1].Tmax*satw ) AC->Tcmd[1] = -AC->Whl[1].Tmax*satw;
        if ( AC->Tcmd[2] >  AC->Whl[2].Tmax*satw ) AC->Tcmd[2] =  AC->Whl[2].Tmax*satw;
        if ( AC->Tcmd[2] < -AC->Whl[2].Tmax*satw ) AC->Tcmd[2] = -AC->Whl[2].Tmax*satw;
  
     }
  
     //retval = 0;
  
     return retval;
  
  }


  int adcsRwTriadTLEh(struct AcType *AC)
  {

     int sunwheel = 0;
     int nonsmooth = 0;

     static int secondsmis = 0, secondsnmis = 0;
     static double chatx = 0., chaty = 0., chatz = 0.;

     int secondsmistop = 1200;
     int secondsnmistop = 24000;
     int transitiontop = 120;
  
     secondsnmis++;
     if (secondsnmis>secondsnmistop) {
        secondsmis++;
     }
     if (secondsmis>secondsmistop) {
        secondsmis = 0;
        secondsnmis = 0;
     }
  
     int retval = 0.;
  
     struct AcCfsCtrlType *C;
     double L1[3],L2[3],L3[3];
     double HxB[3];
     double AngErr;
     long i,j;
  
     C = &AC->CfsCtrl;
  
     if (C->Init) {
        C->Init = 0;
        C->Kunl = 1.0E6;
     }
  
  /* .. Commanded Attitude */
     double orbit2sar = 0.*3.1416/180.;
     double sarend=orbit2sar;
  
     static double sarendf=0.*3.1416/180.;
     static double exsarendf=0.*3.1416/180.;
     static double phif=0., exphif=0.;
  
     exsarendf=sarendf;
     if (secondsmis>transitiontop) sarendf = 0.4*sarend + 0.6*sarendf;  // 0.34
     else sarendf = sarend;
     
     // Non-smooth
     if (nonsmooth==1) 
        sarendf = sarend;
     
     double orbitend=0.;
     if (secondsmis>secondsmistop-transitiontop) sarendf = orbit2sar;
  
     static double phi;
     static double phiq=0.;
     double dphiq=phiq;
     int sigmoid = 0;
  
     //static double dphiq=0.;
     if (AC->GPS[0].Valid) {
        static double exqln[4]={0.,0.,0.,1.};
        CopyUnitV(AC->PosN,L3);
        VxV(AC->PosN,AC->VelN,L2);
        UNITV(L2);
        UNITV(L3);
        VxV(L2,L3,L1);
        UNITV(L1);                                   // TNR
        if (secondsmis==0) phi=orbitend;
        else {
           if (secondsmis<=transitiontop) phi+=orbit2sar/transitiontop;
           if (secondsmis>secondsmistop-transitiontop) phi-=orbit2sar/transitiontop;
           if (secondsmis>transitiontop && secondsmis<=secondsmistop-transitiontop) phi=sarendf;
           // Alternative transition
           if (sigmoid==1) {
              if (secondsmis<=transitiontop) {
                 //phiq=(sin((2.*secondsmis-transitiontop)*0.5*3.1416/transitiontop)+1)*0.5*(sarend-orbitend)+orbitend;
                 dphiq=1.2733*(1+cos((2.*secondsmis-transitiontop)*0.5*2*3.1416/transitiontop))*(sarend-orbitend)/transitiontop*3.1416/2.;
                 phiq+=dphiq*AC->DT;
              } else if (secondsmis>=secondsmistop-transitiontop) {
                 phiq=(sin((2.*(secondsmis-(secondsmistop-transitiontop))-transitiontop)*0.5*3.1416/transitiontop)+1)*0.5*(orbitend-sarend)+sarend;
                 dphiq=cos((2.*(secondsmis-(secondsmistop-transitiontop))-transitiontop)*0.5*3.1416/transitiontop)*(orbitend-sarend)/transitiontop*3.1416/2.;
              } else {
                 phiq=sarendf;
                 dphiq=0.;
              }
           }
        }
  
        // Filter on phi
        exphif = phif;
        if (secondsmis>0) phif = -0.25*phi + 0.75*phif; // -0.34
        else phif = -phi;
  
        // Non-smooth guidance
        if (nonsmooth==1) {
           phif = phi;
           exphif = phi;
        }
  
        for(i=0;i<3;i++) {
           AC->CLN[0][i] = L1[i];
           AC->CLN[1][i] = L2[i]*cos(phif)-L3[i]*sin(phif);
           AC->CLN[2][i] = L3[i]*cos(phif)+L2[i]*sin(phif);
           if (sigmoid==1) {
              AC->CLN[1][i] = L2[i]*cos(phiq)-L3[i]*sin(phiq);
              AC->CLN[2][i] = L3[i]*cos(phiq)+L2[i]*sin(phiq);
           }
        }
        C2Q(AC->CLN,AC->qln);
        AC->wln[0] = 0.;
        AC->wln[1] = MAGV(AC->VelN)/MAGV(AC->PosN)*cos(phif);
        AC->wln[2] = MAGV(AC->VelN)/MAGV(AC->PosN)*sin(phif);
        if (sigmoid==0) {
           if (secondsmis<=transitiontop && secondsnmis>=secondsnmistop) AC->wln[0] -= 2.*orbit2sar/transitiontop;
           else if (secondsmis>=secondsmistop-transitiontop) AC->wln[0] += 2.*orbit2sar/transitiontop;
           else AC->wln[0] += (phif-exphif)/AC->DT;
        } else {
           // Alternative transition
           AC->wln[0] += dphiq; 
        }
        if (AC->qln[0]*exqln[0]+AC->qln[1]*exqln[1]+AC->qln[2]*exqln[2]+AC->qln[3]*exqln[3]<-0.99) {
           AC->qln[0]=-AC->qln[0];
           AC->qln[1]=-AC->qln[1];
           AC->qln[2]=-AC->qln[2];
           AC->qln[3]=-AC->qln[3];
        }
        exqln[0] = AC->qln[0];
        exqln[1] = AC->qln[1];
        exqln[2] = AC->qln[2];
        exqln[3] = AC->qln[3];
     } else {
        for(i=0;i<3;i++) {
           for(j=0;j<3;j++) {
              AC->CLN[i][j] = 0.0;
           }
           AC->CLN[i][i] = 1.0;
           AC->qln[i] = 0.0;
           AC->wln[i] = 0.0;
        }
        AC->qln[3] = 1.0;
        printf("\n******************\n GPS INVALID !!! \n********************\n");
     }
  
     /* .. Attitude Control Reference */
     QxQT(AC->qbn,AC->qln,AC->qbr);
  
     if (AC->qbn[0]*AC->qbn[0]+AC->qbn[1]*AC->qbn[1]+AC->qbn[2]*AC->qbn[2]==0) {
        AC->qbr[0]=0.;
        AC->qbr[1]=0.;
        AC->qbr[2]=0.;
        AC->qbr[3]=1.;
     }
     static double exqbr[4] = {0.,0.,0.,1.};
     if (AC->qbr[0]*exqbr[0]+AC->qbr[1]*exqbr[1]+AC->qbr[2]*exqbr[2]+AC->qbr[3]*exqbr[3]<-0.99) {
        AC->qbr[0]=-AC->qbr[0];
        AC->qbr[1]=-AC->qbr[1];
        AC->qbr[2]=-AC->qbr[2];
        AC->qbr[3]=-AC->qbr[3];
     }
     exqbr[0] = AC->qbr[0];
     exqbr[1] = AC->qbr[1];
     exqbr[2] = AC->qbr[2];
     exqbr[3] = AC->qbr[3];
  
  
     double mbvb = MAGV(AC->bvb);        // Norm of B
     int coriolis = 0;                   // == 1 => coriolis compensation
     double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};   // Inertia Matrix Diagonal
  
     double qc[4] = {AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3]};
     RECTIFYQ(qc);
  
     static double sqc[4]  = {0., 0., 0., 1.};
  
     //if (mbvb > 0.) 
     {
        double m[3], u[3];
        double ki[3] = {0.,0.,0.};
        for (int i = 0; i < 3; i++) {
           if (secondsmis>transitiontop*1. && secondsmis<secondsmistop-transitiontop) sqc[i] = 0.05*qc[i] + 0.95*sqc[i];
           else sqc[i] = 0.;
        }
        for(i=0;i<3;i++)
           FindPDGains(AC->MOI[i][i],0.035/2,0.9, &C->Kr[i],&C->Kp[i]);
        if (secondsmis>0)
           for(i=0;i<3;i++)
              FindPDGains(AC->MOI[i][i],0.35,0.9, &C->Kr[i],&C->Kp[i]);
        if (secondsmis>transitiontop && secondsmis<secondsmistop-transitiontop) 
           for(i=0;i<3;i++) {
              FindPDGains(AC->MOI[i][i],0.35/0.5,0.8, &C->Kr[i],&C->Kp[i]);
              if (i==0) FindPDGains(AC->MOI[i][i],0.35/0.7,0.8, &C->Kr[i],&C->Kp[i]);
           }
        if (secondsnmis<600)
           for(i=0;i<3;i++)
              FindPDGains(AC->MOI[i][i],0.02,1.0, &C->Kr[i],&C->Kp[i]);
        for(i=0;i<3;i++) {
           if (secondsmis>transitiontop*1.0 && secondsmis<secondsmistop-transitiontop) {ki[0]=0.;ki[1]=4.;ki[2]=6.;}
           AC->Tcmd[i] = -C->Kr[i]*AC->wbn[i]-C->Kp[i]*(2.0*(qc[i]+sqc[i]*ki[i]));
        }
        printf("\n ** 2*C->Kp[1] = %f",2*C->Kp[1]);
        /*
        for (int i=0;i<3;i++) {
           double rho = 1.;
           while (fabs(AC->Tcmd[i]) > AC->Whl[i].Tmax && rho>0.05) {
              rho = rho*0.97;
              AC->Tcmd[i] =  -rho*C->Kr[i]*AC->wbn[i]-C->Kp[i]*(2.0*(qc[i]+sqc[i]*ki[i]));
           }      
        }
        */
        double rho = 1.;
        while ((fabs(AC->Tcmd[0]) > AC->Whl[0].Tmax || fabs(AC->Tcmd[1]) > AC->Whl[1].Tmax || fabs(AC->Tcmd[2]) > AC->Whl[2].Tmax) && rho>0.05) {
           rho = rho*0.97;
           for (int i=0;i<3;i++) AC->Tcmd[i] =  -rho*( C->Kr[i]*AC->wbn[i] + C->Kp[i]*(2.0*(qc[i]+sqc[i]*ki[i])) );
        }      
  
        /* .. Momentum Management */
        C->Kunl = 100;
        double um[3]={0.,0.,0.}, mm[3]={0.,0.,0.}, tm[3]={0.,0.,0.};
  
        if (mbvb > 0.) {
           um[0] -= C->Kunl*AC->Whl[0].H;
           um[1] -= C->Kunl*AC->Whl[1].H;
           um[2] -= C->Kunl*AC->Whl[2].H;
           mm[0] = -(um[1]*AC->bvb[2]-um[2]*AC->bvb[1]);
           mm[1] = -(um[2]*AC->bvb[0]-um[0]*AC->bvb[2]);
           mm[2] = -(um[0]*AC->bvb[1]-um[1]*AC->bvb[0]);
           mm[0]/=(2000000.*mbvb*mbvb);mm[1]/=(2000000.*mbvb*mbvb);mm[2]/=(2000000.*mbvb*mbvb);
           tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
           tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
           tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
           AC->Tcmd[0] -= tm[0];
           AC->Tcmd[1] -= tm[1];
           AC->Tcmd[2] -= tm[2];
           if (coriolis==1) {
              double cor[3] = { (J[1]-J[2])*AC->wbn[1]*AC->wbn[2], (J[2]-J[0])*AC->wbn[2]*AC->wbn[0], (J[0]-J[1])*AC->wbn[0]*AC->wbn[1]};
              AC->Tcmd[0] -= cor[0];
              AC->Tcmd[1] -= cor[1];
              AC->Tcmd[2] -= cor[2];
           }
           if (mm[0]> AC->MTB[0].Mmax) mm[0] = AC->MTB[0].Mmax; 
           if (mm[0]<-AC->MTB[0].Mmax) mm[0] = -AC->MTB[0].Mmax;
           if (mm[1]> AC->MTB[1].Mmax) mm[1] = AC->MTB[1].Mmax; 
           if (mm[1]<-AC->MTB[1].Mmax) mm[1] = -AC->MTB[1].Mmax;
           if (mm[2]> AC->MTB[2].Mmax) mm[2] = AC->MTB[2].Mmax; 
           if (mm[2]<-AC->MTB[2].Mmax) mm[2] = -AC->MTB[2].Mmax; 
        }
  
        for (int i = 0; i < 3; i++) {
           AC->Mcmd[i] = mm[i];
        }

        // Adaptation of drag model
        double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
        static int init=1;
        if (init==1) {
           init = 0;
           chatx = chat0;
           chaty = chat0;
           chatz = chat0;
        }
        if (init==0) {
           double pchat = 1E-3;
           chatx = chatx + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[0]*AC->wbn[0];
           chaty = chaty + pchat*AC->Whl[1].H*fabs(AC->Whl[1].H)*J[1]*AC->wbn[1];
           chatz = chatz + pchat*AC->Whl[2].H*fabs(AC->Whl[2].H)*J[2]*AC->wbn[2];
           if (chatx>1.15*chat0) chatx=1.15*chat0;
           if (chatx<0.85*chat0) chatx=0.85*chat0;
           if (chaty>1.15*chat0) chaty=1.15*chat0;
           if (chaty<0.85*chat0) chaty=0.85*chat0;
           if (chatz>1.15*chat0) chatz=1.15*chat0;
           if (chatz<0.85*chat0) chatz=0.85*chat0;
           AC->Tcmd[0] -= chatx*AC->Whl[0].H*fabs(AC->Whl[0].H);
           AC->Tcmd[1] -= chaty*AC->Whl[1].H*fabs(AC->Whl[1].H);
           AC->Tcmd[2] -= chatz*AC->Whl[2].H*fabs(AC->Whl[2].H);
        }

        // Saturation of wheel
        double satw = 1.0;
        if ( AC->Tcmd[0] >  AC->Whl[0].Tmax*satw ) AC->Tcmd[0] =  AC->Whl[0].Tmax*satw;
        if ( AC->Tcmd[0] < -AC->Whl[0].Tmax*satw ) AC->Tcmd[0] = -AC->Whl[0].Tmax*satw;
        if ( AC->Tcmd[1] >  AC->Whl[1].Tmax*satw ) AC->Tcmd[1] =  AC->Whl[1].Tmax*satw;
        if ( AC->Tcmd[1] < -AC->Whl[1].Tmax*satw ) AC->Tcmd[1] = -AC->Whl[1].Tmax*satw;
        if ( AC->Tcmd[2] >  AC->Whl[2].Tmax*satw ) AC->Tcmd[2] =  AC->Whl[2].Tmax*satw;
        if ( AC->Tcmd[2] < -AC->Whl[2].Tmax*satw ) AC->Tcmd[2] = -AC->Whl[2].Tmax*satw;

        //printf("\n -- ORBIT -- um[%d] = %f    AC->Tcmd[%d] = %f ", 0, um[0], 0, AC->Tcmd[0]);
        //printf("\n          -- um[%d] = %f    AC->Tcmd[%d] = %f ", 1, um[1], 1, AC->Tcmd[1]);
        //printf("\n          -- um[%d] = %f    AC->Tcmd[%d] = %f ", 2, um[2], 2, AC->Tcmd[2]);

     }

     static int first=1;
     FILE *FilePtr;
     //if (secondsmis>0) 
     { 
        if (first) {
              first=0;
              FilePtr = fopen("TSAT/missionh.m", "w");
              fprintf(FilePtr, "vm=   [%f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[0], AC->Tcmd[1], AC->Tcmd[2], phi, secondsmis, secondsnmis, phiq, dphiq, AC->qln[0], AC->qln[1], AC->qln[2], AC->qln[3], AC->qbn[0], AC->qbn[1], AC->qbn[2], AC->qbn[3], chatx, chaty, chatz, sarendf, phif, AC->wln[0], sarend);
              fclose(FilePtr);
        } else {
              FilePtr = fopen("TSAT/missionh.m", "a");
              fprintf(FilePtr, "vm=[vm;[%f %f %f %f %f %f %f %f %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f]];\n", AC->qbr[0], AC->qbr[1], AC->qbr[2], AC->qbr[3], AC->Tcmd[0], AC->Tcmd[1], AC->Tcmd[2], phi, secondsmis, secondsnmis, phiq, dphiq, AC->qln[0], AC->qln[1], AC->qln[2], AC->qln[3], AC->qbn[0], AC->qbn[1], AC->qbn[2], AC->qbn[3], chatx, chaty, chatz, sarendf, phif, AC->wln[0], sarend);
              fclose(FilePtr);
        }
     }
  
     //printf("\n k = ( %f, %f, %f, %f, %f, %f )",a_0, e_0, i_0*180/pi, W_0*180/pi, w_0*180/pi, t_0*180/pi);
     //printf("\n e_0 = %f", e_0);
     //printf("\n W_0 = %f", W_0*180/pi);
  
     return retval;
  
  }
  