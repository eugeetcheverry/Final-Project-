#include "42.h"
#include <math.h>

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
   static int secondspointed = 0, secondssun = 0;
   
   double mbvb = MAGV(AC->bvb);        // Norm of B
   static double kq = 0., kw = 10.;    // Initial Gains
   int coriolis = 0;                   // == 1 => coriolis compensation
   int iwearth = 1;                   // >= 0 => M3, axis of rotation, x=0 y=1 z=2
   static double earthvel = 0.;        // != 0 => M3,                   earthvel * wo
   double kw0 = 20.0;                  // Derivative Gain with Eclipse
   double kw1 = 5.0;                   // Derivative Gain without Eclipse
   double kq0 = 0.075;                 // Proportional Gain with Eclipse
   double kq1 = 0.0075;                // Proportional Gain without Eclipse
   double eps = 0.001;                 // Averaging Parameter
   double thnw = 0.001;                // End of Detumbling Velocity Condition [rad/sec] 
   double nm=0.;                       // Auxiliary Variable
   double threclipse = 1.0;            // Eclipse detection threshold (between 0 and 6, 6 never detects daylight, 0 never detects eclipse)
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
         /* Sun Pointing Feedback */
         double sumillum = AC->CSS[0].Illum+AC->CSS[1].Illum+AC->CSS[2].Illum+AC->CSS[3].Illum+AC->CSS[4].Illum+AC->CSS[5].Illum;
         if (sumillum>threclipse) {
            /* Sun with AC->svb */
            secondssun++;
            double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
            double qc[4] = {-svb[2], 0., svb[0], svb[1]};
            static int secondschange=0;
            if (svb[1]<0. && secondschange>100*60) {
               qc[0] = -qc[0];
               qc[1] = -qc[1];
               qc[2] = -qc[2];
               secondschange=0;
            } else {
               secondschange++;
            }
            if (nw < thnw && kq==0.) {kq = kq0; kw = kw0;}
            if (kq!=0.) {
               if (secondssun>100*60) {
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
      /* m = - b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      /* Direct Control Allocation */
      if (fabs(m[0])>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/nm;
         m[1]= m[1]*AC->MTB[0].Mmax/nm;
         m[2]= m[2]*AC->MTB[0].Mmax/nm;
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[1])>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/nm;
         m[1]= m[1]*AC->MTB[1].Mmax/nm;
         m[2]= m[2]*AC->MTB[1].Mmax/nm;
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[2])>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/nm;
         m[1]= m[1]*AC->MTB[2].Mmax/nm;
         m[2]= m[2]*AC->MTB[2].Mmax/nm;
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
      if (secondssun<3000 && secondspointed>500) {
         retval = 1;
      }
      
      // Adaptation of drag model
      if (iwheel>=0 && kq!=0.) {
         double chat0 = AC->Whl[0].Tmax/AC->Whl[0].Hmax/AC->Whl[0].Hmax;
         static double chat = 0.00023/0.00177/0.00177;
         double pchat = 1E8;
         chat = chat + pchat*AC->Whl[0].H*fabs(AC->Whl[0].H)*J[iwheel]*AC->wbn[iwheel];
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
   /*if ( secondssun > 100000) {
      earthvel = -1.0;
   }*/

   //retval = 0.;

   return retval;

}


int adcsRwTriadTLE(struct AcType *AC)
{

   int sunwheel = 0;

   int retval = 0.;
   static int secondspointed = 0, secondssun = 0;

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
      if (AC->GPS[0].Valid) {
         CopyUnitV(AC->PosN,L3);
         VxV(AC->PosN,AC->VelN,L2);
         UNITV(L2);
         UNITV(L3);
         for(i=0;i<3;i++) {
            L2[i] = -L2[i];
            L3[i] = -L3[i];
         }
         VxV(L2,L3,L1);
         UNITV(L1);
         for(i=0;i<3;i++) {
            AC->CLN[0][i] = L1[i];
            AC->CLN[1][i] = L2[i];
            AC->CLN[2][i] = L3[i];
         }
         C2Q(AC->CLN,AC->qln);
         AC->wln[1] = -MAGV(AC->VelN)/MAGV(AC->PosN);
      }
      else {
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
            
/* .. Attitude Control */
      if (AC->StValid) {
         QxQT(AC->qbn,AC->qln,AC->qbr);
         RECTIFYQ(AC->qbr);
      }
      else {
         for(i=0;i<3;i++) AC->qbr[i] = 0.0;
         AC->qbr[3] = 1.0;
      }
      double kw = 30.0, kq = 0.15, eps=0.001;
      double svb[3] = {AC->svb[0], AC->svb[1], AC->svb[2]};
      double nsvb = sqrt(svb[0]*svb[0]+svb[1]*svb[1]+svb[2]*svb[2]);
      double qc[4] = {-svb[2]/nsvb, 0., svb[0]/nsvb, 0.};
      double J[3] = {AC->MOI[0][0], AC->MOI[1][1], AC->MOI[2][2]};
      for(i=0;i<3;i++) {
         C->therr[i] = Limit(2.0*AC->qbr[i],-0.1,0.1);
         C->werr[i] = AC->wbn[i] - AC->wln[i];
         AC->Tcmd[i] = Limit(-C->Kr[i]*C->werr[i] - C->Kp[i]*C->therr[i],-0.1,0.1);
         if (sunwheel) {
            C->therr[i] = qc[i];
            C->werr[i] = AC->wbn[i];
            AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - kq * eps * eps * (1/J[i]) * qc[i];
         }
      }
/* .. Momentum Management */
      for(i=0;i<3;i++) {
         AC->Hvb[i] = AC->MOI[i][i]*AC->wbn[i];
         for(j=0;j<AC->Nwhl;j++) AC->Hvb[i] += AC->Whl[j].Axis[i]*AC->Whl[j].H;
      }
      VxV(AC->Hvb,AC->bvb,HxB);
      for(i=0;i<3;i++) AC->Mcmd[i] = C->Kunl*HxB[i];

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
         double sumillum = AC->CSS[0].Illum+AC->CSS[1].Illum+AC->CSS[2].Illum+AC->CSS[3].Illum+AC->CSS[4].Illum+AC->CSS[5].Illum;
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
      /* m = - b x u */
      m[0] = -(u[1]*AC->bvb[2]-u[2]*AC->bvb[1]);
      m[1] = -(u[2]*AC->bvb[0]-u[0]*AC->bvb[2]);
      m[2] = -(u[0]*AC->bvb[1]-u[1]*AC->bvb[0]);
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      /* Direct Control Allocation */
      if (fabs(m[0])>AC->MTB[0].Mmax) {
         m[0]= m[0]*AC->MTB[0].Mmax/nm;
         m[1]= m[1]*AC->MTB[0].Mmax/nm;
         m[2]= m[2]*AC->MTB[0].Mmax/nm;
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[1])>AC->MTB[1].Mmax) {
         m[0]= m[0]*AC->MTB[1].Mmax/nm;
         m[1]= m[1]*AC->MTB[1].Mmax/nm;
         m[2]= m[2]*AC->MTB[1].Mmax/nm;
      }
      nm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
      if (fabs(m[2])>AC->MTB[2].Mmax) {
         m[0]= m[0]*AC->MTB[2].Mmax/nm;
         m[1]= m[1]*AC->MTB[2].Mmax/nm;
         m[2]= m[2]*AC->MTB[2].Mmax/nm;
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


         /*for(int i=0;i<3;i++) {
            AC->Hvb[i] = AC->MOI[i][i]*AC->wbn[i];
         }
         VxV(AC->Hvb,AC->bvb,HxB);
         um[0] = -C->Kunl*HxB[0];
         um[1] = -C->Kunl*HxB[1];
         um[2] = -C->Kunl*HxB[2];
         */
            
            /*for mm tm[0] = (mm[1]*AC->bvb[2]-mm[2]*AC->bvb[1]);
            tm[1] = (mm[2]*AC->bvb[0]-mm[0]*AC->bvb[2]);
            tm[2] = (um[0]*AC->bvb[1]-mm[1]*AC->bvb[0]);
            AC->Tcmd[iwheel] -= tm[iwheel];*/


            /*
            double nqp[3] = {-svb[2], 0., svb[0]};
            if (svb[1]<0. && secondssun<100*60) {
               nqp[0]=-nqp[0];
               nqp[1]=-nqp[1];
               nqp[2]=-nqp[2];
            }
            double auxnorm = sqrt(nqp[0]*nqp[0]+nqp[1]*nqp[1]+nqp[2]*nqp[2]);
            double angqp = acos(svb[0]);
            if (auxnorm>0.) {
               qc[0]=nqp[0]/auxnorm*sin(angqp/2);
               qc[1]=nqp[1]/auxnorm*sin(angqp/2);
               qc[2]=nqp[2]/auxnorm*sin(angqp/2);
               qc[3]=cos(angqp/2);
            } else {
               qc[0]=0.;
               qc[1]=0.;
               qc[2]=0.;
               qc[3]=1.;
               if (svb[1]<0.) qc[3]=-1.;
            }
            */

/*

            double lvlh=0.;
            if (lvlh==1.) {
               double L1[3],L2[3],L3[3];
               if (AC->GPS[0].Valid) {
                  CopyUnitV(AC->PosN,L3);
                  VxV(AC->PosN,AC->VelN,L2);
                  UNITV(L2);
                  UNITV(L3);
                  for(i=0;i<3;i++) {
                     L2[i] = -L2[i];
                     L3[i] = -L3[i];
                  }
                  VxV(L2,L3,L1);
                  UNITV(L1);
                  for(i=0;i<3;i++) {
                     AC->CLN[0][i] = L1[i];
                     AC->CLN[1][i] = L2[i];
                     AC->CLN[2][i] = L3[i];
                  }
                  C2Q(AC->CLN,AC->qln);
                  AC->wln[1] = -MAGV(AC->VelN)/MAGV(AC->PosN);
               } else {
                  for(i=0;i<3;i++) {
                     for(int j=0;j<3;j++) {
                        AC->CLN[i][j] = 0.0;
                     }
                     AC->CLN[i][i] = 1.0;
                     AC->qln[i] = 0.0;
                     AC->wln[i] = 0.0;
                  }
                  AC->qln[3] = 1.0;
               }
               if (AC->StValid) {
                  QxQT(AC->qbn,AC->qln,AC->qbr);
                  RECTIFYQ(AC->qbr);
               } else {
                  for(i=0;i<3;i++) AC->qbr[i] = 0.0;
                  AC->qbr[3] = 1.0;
               }
               qc[0]=AC->qbr[0];
               qc[1]=AC->qbr[1];
               qc[2]=AC->qbr[2];
            }


*/

/*

            double nq = sqrt(svb[0]*svb[0]+svb[2]*svb[2]);
            if (nq==0.) nq=0.00001;
            double angleq = asin(nq);
            double ss = 1;if (svb[1]<0.) ss = -1;
            double qc[4] = {-svb[2]/nq*sin(angleq/2.), 0., svb[0]/nq*sin(angleq/2.), ss*cos(angleq/2.)};

*/

/*
            double nq = sqrt(svb[0]*svb[0]+svb[2]*svb[2]);
            double angleq = asin(nq);
            double ss = 1;if (svb[1]<0.) ss = -1;
            if (nq==0.) nq=1.;
            double qc[4] = {-svb[2]/nq*sin(angleq/2.), 0., svb[0]/nq*sin(angleq/2.), ss*cos(angleq/2.)};
            // True quaternion compensated *2
            qc[0]=2*qc[0];qc[2]=2*qc[2];
            // Simple quaternion
            qc[0]=-svb[2];qc[2]=svb[0];

*/

/*

         static int first=1;
         FILE *FilePtr;
         if (first) {
            first=0;
            FilePtr = fopen("TSAT/chat.m", "w");
            fprintf(FilePtr, "vc=[%f];\n", chat);
            fclose(FilePtr);
         } else {
            FilePtr = fopen("TSAT/chat.m", "a");
            fprintf(FilePtr, "vc=[vc;%f];\n", chat);
            fclose(FilePtr);
         }


*/

/*

   // Earth Pointing Velocity Test
   // To be sent by command, this is just an example
   if ( secondssun>300000 ) {
      if (AC->svb[1]>0.) {
         earthvel=0.25; 
      } else { 
         earthvel=-0.25;
      }
   } 
   if ( secondssun>500000 ) {
      if (AC->svb[1]>0.) {
         earthvel=0.5; 
      } else { 
         earthvel=-0.5;
      }
   }

*/


/*

   static int mode4 = 0;
   double rollref=0.;

*/

/*
            // Mode 4
            if (mode4!=0 && secondssun>5*60) {
               // on = s x b
               double on[3];
               on[0] = -(svb[1]*AC->bvb[2]-svb[2]*AC->bvb[1]);
               on[1] = -(svb[2]*AC->bvb[0]-svb[0]*AC->bvb[2]);
               on[2] = -(svb[0]*AC->bvb[1]-svb[1]*AC->bvb[0]);
               double non = sqrt(on[0]*on[0]+on[1]*on[1]+on[2]*on[2]);
               on[0] = on[0]/non;
               on[1] = on[1]/non;
               on[2] = on[2]/non;
               if (mode4==1) rollref = 0.01*on[0];
               if (mode4==2) rollref = -0.1*on[2];
               if (svb[1]<0.) rollref = -rollref;
               printf("rollref=%f   ssun=%d  ", rollref, secondssun);
            }
            // end of Mode 4
*/

/*

After AC->Tcmd[i] = - kw * eps * J[i] * AC->wbn[i] - kq * eps * eps * (1/J[i]) * qc[i];

               if (mode4!=0) AC->Tcmd[i] = AC->Tcmd[i] - kq * eps * eps * (1/J[i]) * rollref;

*/



   /*
   // Mode 4 Test
   // To be sent by command, this is just an example for noon
   if ( secondssun>5*60 ) {
      mode4 = 2;
   }
   static int first=1;
   FILE *FilePtr;
   if (first) {
            first=0;
            FilePtr = fopen("TSAT/rref.m", "w");
            fprintf(FilePtr, "vr=[%f];\n", asin(rollref*9.98)*180./3.1416);
            fclose(FilePtr);
   } else {
            FilePtr = fopen("TSAT/rref.m", "a");
            fprintf(FilePtr, "vr=[vr;%f];\n", asin(rollref*9.98)*180./3.1416);
            fclose(FilePtr);
   }
   */
