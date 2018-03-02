#include "ChargedParticle.h"
#include <cmath>
#include <iostream>
#include "mymath.h"
#include "EFieldMap.h"
#include "EFieldPoint.h"
#include "MagFieldMap.h"
#include "Matrix.h"
#include "globals.h"

ChargedParticle::ChargedParticle()
{
  DefineConstants();
}

ChargedParticle::ChargedParticle(const Vector3 & pos,
				 const Vector3 & vel,
				 const double & cha,
				 const double & ma,
				 const double & potent,
				 const double & ti,
				 const double & sre,
				 const double & srb)
{
  position = pos;
  oldposition = pos;
  velocity = vel;
  charge = cha;
  mass = ma;
  pot_E = potent;
  kin_E = 0.5 * mass * velocity.dot(velocity);
  time = ti;

  oldvelocity = vel;
  oldB = NullVector;

  c2m_r = charge / mass;
  
  search_radiuse = sre;
  search_radiusb = srb;
  DefineConstants();
}

ChargedParticle::~ChargedParticle()
{

}

void ChargedParticle::DefineConstants()
{
  c_rkf45[0]=0;
  c_rkf45[1]=0;
  c_rkf45[2]=1.0/4.0;
  c_rkf45[3]=3.0/8.0;
  c_rkf45[4]=12.0/13.0;
  c_rkf45[5]=1.0;
  c_rkf45[6]=0.5;
 
  c_rkkc45[0]=0;
  c_rkkc45[1]=0;
  c_rkkc45[2]=1.0/5.0;
  c_rkkc45[3]=3.0/10.0;
  c_rkkc45[4]=3.0/5.0;
  c_rkkc45[5]=1.0;
  c_rkkc45[6]=7.0/8.0;
  
  for(int i= 0; i < 7; i++) {
    for(int j= 0; j < 7; j++) {
      a_rkf45[i][j] = 0;
      a_rkkc45[i][j] = 0;
    }
    b_rkf451[i] = 0;
    b_rkf452[i] = 0;
    b_rkkc451[i] = 0;
    b_rkkc452[i] = 0;
  }
  
  a_rkf45[2][1] = 1.0/4.0;
  a_rkf45[3][1] = 3.0/32.0;
  a_rkf45[3][2] = 9.0/32.0;
  a_rkf45[4][1] = 1932.0/2197.0;
  a_rkf45[4][2] = -7200.0/2197.0;
  a_rkf45[4][3] = 7296.0/2197.0;
  a_rkf45[5][1] = 439.0/216.0;
  a_rkf45[5][2] = -8.0;
  a_rkf45[5][3] = 3680.0/513.0;
  a_rkf45[5][4] = -845.0/4104.0;
  a_rkf45[6][1] = -8.0/27.0;
  a_rkf45[6][2] = 2.0;
  a_rkf45[6][3] = -3544.0/2565.0;
  a_rkf45[6][4] = 1859.0/4104.0;
  a_rkf45[6][5] = -11.0/40.0;

  b_rkf451[0]=0;
  b_rkf451[1]=25.0/216;
  b_rkf451[2]=0;
  b_rkf451[3]=1408.0/2565.0;
  b_rkf451[4]=2197.0/4104.0;
  b_rkf451[5]=-1.0/5.0;
  b_rkf451[6]=0;

  b_rkf452[0]=0;
  b_rkf452[1]=16.0/135.0;
  b_rkf452[2]=0;
  b_rkf452[3]=6656.0/12825.0;
  b_rkf452[4]=28561.0/56430.0;
  b_rkf452[5]=-9.0/50.0;
  b_rkf452[6]=2.0/55.0;
  
  a_rkkc45[2][1] = 1.0/5.0;
  a_rkkc45[3][1] = 3.0/40.0;
  a_rkkc45[3][2] = 9.0/40.0;
  a_rkkc45[4][1] = 3.0/10.0;
  a_rkkc45[4][2] = -9.0/10.0;
  a_rkkc45[4][3] = 6.0/5.0;
  a_rkkc45[5][1] = -11.0/54.0;
  a_rkkc45[5][2] = 5.0/2.0;
  a_rkkc45[5][3] = -70.0/27.0;
  a_rkkc45[5][4] = 35.0/27.0;
  a_rkkc45[6][1] = 1631.0/55296.0;
  a_rkkc45[6][2] = 175.0/512.0;
  a_rkkc45[6][3] = 575.0/13828.0;
  a_rkkc45[6][4] = 44275.0/110592.0;
  a_rkkc45[6][5] = 253.0/4096.0;

  b_rkkc451[0]=0;
  b_rkkc451[1]=37.0/378.0;
  b_rkkc451[2]=0;
  b_rkkc451[3]=250.0/621.0;
  b_rkkc451[4]=125.0/594.0;
  b_rkkc451[5]=0;
  b_rkkc451[6]=512.0/1771.0;

  b_rkkc452[0]=0;
  b_rkkc452[1]=2825.0/27648.0;
  b_rkkc452[2]=0;
  b_rkkc452[3]=18575.0/48384.0;
  b_rkkc452[4]=13525.0/55296.0;
  b_rkkc452[5]=277.0/14336.0;
  b_rkkc452[6]=1.0/4.0;

  c_dp[0]=0;
  c_dp[1]=0;
  c_dp[2]=1/18;
  c_dp[3]=1/12;
  c_dp[4]=1/8;
  c_dp[5]=5/16;
  c_dp[6]=3/8;
  c_dp[7]=59/400;
  c_dp[8]=93/200;
  c_dp[9]=5490023248/9719169821;
  c_dp[10]=13/20;
  c_dp[11]=30992876149296355/33518267164510641;
  c_dp[12]=1;
  c_dp[13]=1;
  
  for(int i= 0; i < 14; i++) {
    for(int j= 0; j < 14; j++) {
      a_dp[i][j] = 0;
    }
    b_dp1[i] = 0;
    b_dp2[i] = 0;
  }

  a_dp[2][1]=.5555555555555555555555555555555555555555555555555555555555555555555555555555555555556e-1;
  a_dp[3][1]=.2083333333333333333333333333333333333333333333333333333333333333333333333333333333333e-1;
  a_dp[3][2]=.625e-1;
  a_dp[4][1]=.3125e-1;
  a_dp[4][2]=0.;
  a_dp[4][3]=.9375e-1;
  a_dp[5][1]=.3125;
  a_dp[5][2]=0.;
  a_dp[5][3]=-1.171875;
  a_dp[5][4]=1.171875;
  a_dp[6][1]=.375e-1;
  a_dp[6][2]=0.;
  a_dp[6][3]=0.;
  a_dp[6][4]=.1875;
  a_dp[6][5]=.15;
  a_dp[7][1]=.4791013711111111111111111111111111111111111111111111111111111111111111111111111111111e-1;
  a_dp[7][2]=0.;
  a_dp[7][3]=0.;
  a_dp[7][4]=.1122487127777777777777777777777777777777777777777777777777777777777777777777777777778;
  a_dp[7][5]=-.2550567377777777777777777777777777777777777777777777777777777777777777777777777777778e-1;
  a_dp[7][6]=.1284682388888888888888888888888888888888888888888888888888888888888888888888888888889e-1;
  a_dp[8][1]=.1691798978729228118143110713603823606551492879543441068765183916901985069878116840258e-1;
  a_dp[8][2]=0.;
  a_dp[8][3]=0.;
  a_dp[8][4]=.3878482784860431695265457441593733533707275558526027137301504136188413357699752956243;
  a_dp[8][5]=.3597736985150032789670088963477236800815873945874968299566918115320174598617642813621e-1;
  a_dp[8][6]=.1969702142156660601567152560721498881281698021684239074332818272245101975612112344987;
  a_dp[8][7]=-.1727138523405018387613929970023338455725710262752107148467532611655731300161441266618;
  a_dp[9][1]=.6909575335919230064856454898454767856104137244685570409618590205329379141801779261271e-1;
  a_dp[9][2]=0.;
  a_dp[9][3]=0.;
  a_dp[9][4]=-.6342479767288541518828078749717385465426336031120747443029278238678259156558727343870;
  a_dp[9][5]=-.1611975752246040803668769239818171234422237864808493434053355718725564004275765826294;
  a_dp[9][6]=.1386503094588252554198669501330158019276654889495806914244800957316809692178564181463;
  a_dp[9][7]=.9409286140357562697242396841302583434711811483145028697758469500622977999086075976730;
  a_dp[9][8]=.2116363264819439818553721171319021047635363886083981439225363020792869599016568720713;
  a_dp[10][1]=.1835569968390453854898060235368803084497642516277329033567822939550366893470341427061;
  a_dp[10][2]=0.;
  a_dp[10][3]=0.;
  a_dp[10][4]=-2.468768084315592452744315759974107457777649753547732495587510208118948846864075671103;
  a_dp[10][5]=-.2912868878163004563880025728039519800543380294081727678722180378521228988619909525814;
  a_dp[10][6]=-.2647302023311737568843979946594614325963055282676866851254669985184857900430792761417e-1;
  a_dp[10][7]=2.847838764192800449164518254216773770231582185011953158945518598604677368771006006947;
  a_dp[10][8]=.2813873314698497925394036418267117820980705455360140173438806414700945017562775967023;
  a_dp[10][9]=.1237448998633146576270302126636397203122013536069738523260934117931117648560568049432;
  a_dp[11][1]=-1.215424817395888059160510525029662994880024988044112261492933435562347094075572196306;
  a_dp[11][2]=0.;
  a_dp[11][3]=0.;
  a_dp[11][4]=16.67260866594577243228041328856410774858907460078758981907659463146817850515953150318;
  a_dp[11][5]=.9157418284168179605957186504507426331593736334298114556675054711691705693180672015549;
  a_dp[11][6]=-6.056605804357470947554505543091634004081967083324413691952297768849489422976536037188;
  a_dp[11][7]=-16.00357359415617811184170641007882303068079304063907122528682690677051264091851070681;
  a_dp[11][8]=14.84930308629766255754539189802663208272299893302745636963476380528935538206408294404;
  a_dp[11][9]=-13.37157573528984931829304139618159579089195289469740214314819172008509426917249413996;
  a_dp[11][10]=5.134182648179637933173253611658602898712494286162881495270691720760048063805245199662;
  a_dp[12][1]=.2588609164382642838157309322317577667296307766301063163257807024364127266506000943372;
  a_dp[12][2]=0.;
  a_dp[12][3]=0.;
  a_dp[12][4]=-4.774485785489205112310117509706042746829391853746564335432052648850273207666153414325;
  a_dp[12][5]=-.4350930137770325094407004118103177819323551661617974116300885410959889587870529265443;
  a_dp[12][6]=-3.049483332072241509560512866312031613982854911220735245095857885009959537455491093975;
  a_dp[12][7]=5.577920039936099117423676634464941858623588944531498679305832747426169795721583100328;
  a_dp[12][8]=6.155831589861040097338689126688954481197754937462937466615262374558053507207577588692;
  a_dp[12][9]=-5.062104586736938370077406433910391644990220712141673880668482097753119682588189892664;
  a_dp[12][10]=2.193926173180679061274914290465806019788262707389033759251111926114017568440058142981;
  a_dp[12][11]=.1346279986593349415357262378873236613955852772571946513284934221746877884770684011715;
  a_dp[13][1]=.8224275996265074779631682047726665909572303617765850630130165537026064642659285212794;
  a_dp[13][2]=0.;
  a_dp[13][3]=0.;
  a_dp[13][4]=-11.65867325727766428397655303545841477547369082638864247596731917545919361630644756876;
  a_dp[13][5]=-.7576221166909361958811161540882449653663757591954118634754438929149012424414834787423;
  a_dp[13][6]=.7139735881595815279782692827650546753142248878566301594716125352788068216039494192152;
  a_dp[13][7]=12.07577498689005673956617044860067967095705800972197897084832370633043082304810801069;
  a_dp[13][8]=-2.127659113920402656390820858969398635427927973275119750336406473203376323696576615278;
  a_dp[13][9]=1.990166207048955418328071698344314152176173017359794982793519250552799420955298804232;
  a_dp[13][10]=-.2342864715440402926602946918568015314512425817004394700255034276765120670064339827205;
  a_dp[13][11]=.1758985777079422650731051058901448183145508638446243836782009233893397195776568900819;
  a_dp[13][12]=0.;


  b_dp1[1]=.4174749114153024622208592846850711513419360280744778485846131359556941056165592859127e-1;
  b_dp1[2]=0.;
  b_dp1[3]=0.;
  b_dp1[4]=0.;
  b_dp1[5]=0.;
  b_dp1[6]=-.5545232861123930896152189465471671889359328156519226303661201143086248409346173220672e-1;
  b_dp1[7]=.2393128072011800970467473542487569696603053593110659071772286813575015789696127063297;
  b_dp1[8]=.7035106694034430230580464108897021513663799729422329988912282232809616439130702896567;
  b_dp1[9]=-.7597596138144609298844876770850584076554230192157088121727371293320966312254783676226;
  b_dp1[10]=.6605630309222863414613785948378206399404197125341442757648140724130090459250911959182;
  b_dp1[11]=.1581874825101233355296148386006854439728567324034642678211136427219891617064005758134;
  b_dp1[12]=-.2381095387528628044718635553056971935251390792174541593034967926060717257568905964800;
  b_dp1[13]=.25;
  b_dp2[1]=.2955321367635349698196488311203224657733277925009198076731622388624965866492863628737e-1;
  b_dp2[2]=0.;
  b_dp2[3]=0.;
  b_dp2[4]=0.;
  b_dp2[5]=0.;
  b_dp2[6]=-.8286062764877970397668056126887191847354037544018122961965538027330129339209559093501;
  b_dp2[7]=.3112409000511183279299137516268570512897363217826855582698888340873138124416908811814;
  b_dp2[8]=2.467345190599886981964685704068761458562259575475834628418174165533324998199993514039;
  b_dp2[9]=-2.546941651841908739127380075415708961780905163114681004322493329883513805107784462867;
  b_dp2[10]=1.443548583676775240301874950690104268511068010189240811834537552274685797398657155566;
  b_dp2[11]=.7941559588112728727130195416222867713146778637419587678468591239050802787902574069868e-1;
  b_dp2[12]=.4444444444444444444444444444444444444444444444444444444444444444444444444444444444444e-1;
  b_dp2[13]=0.;
}

void ChargedParticle::Prop_DT(const Vector3 & theForce, double timestep, const my_kd_tree_t & index)
{
  Vector3 r_0 = position;
  double e_total = kin_E + pot_E;
  double ke_old = 0;
  double pe_old = pot_E;
  Vector3 positionnew = position + velocity * timestep + 
    theForce * (squ(timestep) / (2*mass));
  Vector3 velocitynew = velocity + theForce * (timestep / mass);

  position = positionnew;
  velocity = velocitynew;
  kin_E = 0.5 * mass * velocity.dot(velocity);

}

void ChargedParticle::SetTolerance(double tol)
{
  CTolerance = tol;
}

void ChargedParticle::Prop_DT(const Vector3 & theForce, double timestep,
			      const Vector3 & workForce, int changes, const my_kd_tree_t & index)
{
  Vector3 r_0 = position;
  double e_total = kin_E + pot_E;
  double ke_old = 0;
  double pe_old = pot_E;
  Vector3 positionnew = position + velocity * timestep + 
    theForce * (squ(timestep) / (2*mass));
  Vector3 velocitynew = velocity + theForce * (timestep / mass);

  position = positionnew;

  Vector3 tworkForce;
  double work = tworkForce.dot(r_0 - positionnew);
  if(work < 1e-15) work = 0;
  pot_E = pot_E - work;

  double ke_new = 0;
  // if (velocitynew.xc() != velocity.xc()) {
  //   ke_new += 0.5 * mass * squ(velocitynew.xc());
  //   ke_old += 0.5 * mass * squ(velocity.xc());
  //   changes += 1;
  // }
  // if (velocitynew.yc() != velocity.yc()) {
  //   ke_new += 0.5 * mass * squ(velocitynew.yc());
  //   ke_old += 0.5 * mass * squ(velocity.yc());
  //   changes += 2;
  // }
  // if (velocitynew.zc() != velocity.zc()) {
  //   ke_new += 0.5 * mass * squ(velocitynew.zc());
  //   ke_old += 0.5 * mass * squ(velocity.zc());
  //   changes += 4;
  // }
  if (changes & 1 == 1) {
    ke_new += 0.5 * mass * squ(velocitynew.xc());
    ke_old += 0.5 * mass * squ(velocity.xc());
  }
  if (changes & 2 == 2) {
    ke_new += 0.5 * mass * squ(velocitynew.yc());
    ke_old += 0.5 * mass * squ(velocity.yc());
  }
  if (changes & 4 == 4) {
    ke_new += 0.5 * mass * squ(velocitynew.zc());
    ke_old += 0.5 * mass * squ(velocity.zc());
  }

  double scale = sqrt((ke_old + 2*work/mass)/ke_new);
  double tx, ty, tz;
  velocitynew.return_coords(tx,ty,tz);

  //std::cout << ke_old << ' ' << work << ' ' << ke_new << ' ' << scale << std::endl;
    //std::cout << scale << ' ' << velocity << ' ' << velocitynew << std::endl;

  switch (changes) {
  case 1:
    velocity.set_coords(scale * tx, ty, tz);
    break;
  case 2:
    velocity.set_coords(tx, scale * ty, tz);
    break;
  case 3:
    velocity.set_coords(scale * tx, scale * ty, tz);
    break;
  case 4:
    velocity.set_coords(tx, ty, scale * tz);
    break;
  case 5:
    velocity.set_coords(scale * tx, ty, scale * tz);
    break;
  case 6:
    velocity.set_coords(tx, scale * ty, scale * tz);
    break;
  case 8:
    velocity.set_coords(scale * tx, scale * ty, scale * tz);
    break;
  }

  kin_E = 0.5 * mass * velocity.dot(velocity);

  
  // std::cout << velocity << ' ' << (e_total - pot_E) 
  // 	    << ' ' << kin_E << ' ' << (e_total - pot_E) / kin_E << std::endl;

}

void ChargedParticle::Prop_DT_RK4(double t, const EFieldMap & efield,
				  const MagFieldMap & bfield, const my_kd_tree_t & index)
{
  Vector3 k1[2],k2[2],k3[2],k4[2];
  double h = t;
  k1[0] = h*velocity;
  k1[1] = h*GetForce(position, velocity, t, efield, bfield, index);
  k2[0] = h*(velocity + 0.5*k1[1]);
  k2[1] = h*GetForce(position+k1[0]*0.5,
		     k2[0] / h,          t+h/2.0, efield, bfield, index);
  k3[0] = h*(velocity + 0.5*k2[1]);
  k2[1] = h*GetForce(position+k2[0]*0.5,
		     k3[0] / h,         t+h/2.0, efield, bfield, index);
  k4[0] = h*(velocity + k3[1]);
  k2[1] = h*GetForce(position+k3[0],
		     k4[0] / h,         t+h, efield, bfield, index);
  position = position + (k1[0] + 2.0*(k2[0] + k3[0]) + k4[0])/6.0;
  velocity = velocity + (k1[1] + 2.0*(k2[1] + k3[1]) + k4[1])/6.0;
  time += t;
}

void ChargedParticle::Prop_DT_RKF45(const EFieldMap & efield,
				    const MagFieldMap & bfield, 
				    const my_kd_tree_t & index,
				    double & stepsize)
{
  // Vector3 k1[2],k2[2],k3[2],k4[2],k5[2],k6[2];
  // Vector3 v = velocity;
  // Vector3 x = position;
  // Vector3 v_1,v_2;
  // Vector3 x_1,x_2;
  // double h = stepsize;
  // double dvlen, dxlen;
  // double t = 0;
  // dvlen = dxlen = CTolerance;
  // do {
  //   if(dxlen > CTolerance || dvlen > CTolerance) {
  //     h /=2; //= 0.840896 * sqrt(sqrt(CTolerance * h / (dvlen)));
  //   }
  //   else if(dvlen*10 < CTolerance && dxlen*10 < CTolerance) {
  //     h *=2; //= 0.840896 * sqrt(sqrt(CTolerance * h / (dvlen)));
  //   }
  //   //h = rkf45_stepsize;
  //   k1[0] = h*v;
  //   k1[1] = h*GetForce(x,v,t,efield,bfield, index);
  //   //std::cout << "finished k1[0] k1[1]" << std::endl;

  //   k2[0] = h*(v+k1[1]*a_rkf45[2][1]);
  //   k2[1] = h*GetForce(x+k1[0]*a_rkf45[2][1],
  // 		       k2[0]/h,t+h*c_rkf45[2],efield,bfield, index);

  //   k3[0] = h*(v+k1[1]*a_rkf45[3][1]+
  // 	         k2[1]*a_rkf45[3][2]);
  //   k3[1] = h*GetForce(x+k1[0]*a_rkf45[3][1]+
  // 		         k2[0]*a_rkf45[3][2],
  //   		       k3[0]/h,t+h*c_rkf45[3],efield,bfield, index);

  //   k4[0] = h*(v+k1[1]*a_rkf45[4][1]+
  // 	         k2[1]*a_rkf45[4][2]+
  // 	         k3[1]*a_rkf45[4][3]);
  //   k4[1] = h*GetForce(x+k1[0]*a_rkf45[4][1]+
  // 		         k2[0]*a_rkf45[4][2]+
  // 		         k3[0]*a_rkf45[4][3],
  //   		       k4[0]/h,t+h*c_rkf45[4],efield,bfield, index);

  //   k5[0] = h*(v+k1[1]*a_rkf45[5][1]+
  // 	         k2[1]*a_rkf45[5][2]+
  // 	         k3[1]*a_rkf45[5][3]+
  // 	         k4[1]*a_rkf45[5][4]);
  //   k5[1] = h*GetForce(x+k1[0]*a_rkf45[5][1]+
  // 		         k2[0]*a_rkf45[5][2]+
  // 		         k3[0]*a_rkf45[5][3]+
  // 		         k4[0]*a_rkf45[5][4],
  //   		       k5[0]/h,t+h*c_rkf45[5],efield,bfield, index);

  //   k6[0] = h*(v+k1[1]*a_rkf45[6][1]+
  // 	         k2[1]*a_rkf45[6][2]+
  // 	         k3[1]*a_rkf45[6][3]+
  // 	         k4[1]*a_rkf45[6][4]+
  // 	         k5[1]*a_rkf45[6][5]);
  //   k6[1] = h*GetForce(x+k1[0]*a_rkf45[6][1]+
  // 		         k2[0]*a_rkf45[6][2]+
  // 		         k3[0]*a_rkf45[6][3]+
  // 		         k4[0]*a_rkf45[6][4]+
  // 		         k5[0]*a_rkf45[6][5],
  //   		       k6[0]/h,t+h*c_rkf45[6],efield,bfield, index);

  //   v_1 = v+(k1[1]*b_rkf451[1]+
  // 	     k2[1]*b_rkf451[2]+
  // 	     k3[1]*b_rkf451[3]+
  // 	     k4[1]*b_rkf451[4]+
  // 	     k5[1]*b_rkf451[5]);

  //   v_2 = v+(k1[1]*b_rkf452[1]+
  // 	     k2[1]*b_rkf452[2]+
  // 	     k3[1]*b_rkf452[3]+
  // 	     k4[1]*b_rkf452[4]+
  // 	     k5[1]*b_rkf452[5]+
  // 	     k6[1]*b_rkf452[6]);

  //   x_1 = x+(k1[0]*b_rkf451[1]+
  // 	     k2[0]*b_rkf451[2]+
  // 	     k3[0]*b_rkf451[3]+
  // 	     k4[0]*b_rkf451[4]+
  // 	     k5[0]*b_rkf451[5]);

  //   x_2 = x+(k1[0]*b_rkf452[1]+
  // 	     k2[0]*b_rkf452[2]+
  // 	     k3[0]*b_rkf452[3]+
  // 	     k4[0]*b_rkf452[4]+
  // 	     k5[0]*b_rkf452[5]+
  // 	     k6[0]*b_rkf452[6]);

  //   Vector3 dx = x_1 - x_2;
  //   Vector3 dv = v_1 - v_2;
    
  //   dxlen = sqrt(dx.dot(dx));
  //   dvlen = sqrt(dv.dot(dv))*h;
    
  //   // std::cerr << rkf45_stepsize << ' ' <<
  //   //   dvlen << ' ' << dxlen << std::endl;
  //   //std::cout << "finished all" << std::endl;
  // }while ((dxlen > CTolerance || dvlen > CTolerance));

  // position = x_1;
  // velocity = v_1;

  // time += h;
  //std::cout << "finished all" << std::endl;
}

bool ChargedParticle::Prop_DT_RKF451(const EFieldMap & efield,
				    const EFieldMap & bfield, 
				    const my_kd_tree_t & index_e,
				    const my_kd_tree_t & index_b,
				    double & stepsize)
{
  Vector3 k1[2],k2[2],k3[2],k4[2],k5[2],k6[2];
  Vector3 v = velocity;
  Vector3 x = position;
  Vector3 v_1,v_2;
  Vector3 x_1,x_2;
  double h = stepsize;
  double dvlen, dxlen;
  double t = 0;
  dvlen = dxlen = CTolerance;
  int count = 0;
  
  Vector3 dx = x_1 - x_2;
  Vector3 dv = v_1 - v_2;
  bool good_force = true;

  do {
    count++;
    if(count == 4) {
      h=h;
    }
    else if(dxlen > CTolerance) {
      // std::cerr << h << ' ' << dxlen << ' ' << pow(CTolerance / dxlen, 0.2) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      //  		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2) << " s"
      //  		<< std::endl;
      h = 0.840896 * h * pow(CTolerance / abs(dxlen), 0.2);
      if(count == 4 || h != h) h = 1e-12;
    }
    else if(dxlen < CTolerance/50.0) {
      // std::cerr << h << ' ' << dxlen << ' ' << pow(CTolerance / dxlen, 0.2) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      //  		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2) << " b"
      //  		<< std::endl;
      h = 0.840896 * h * pow(CTolerance / abs(dxlen), 0.25);
      if(count == 4 || h != h || h > 1e-10 || dxlen == 0) h = 1e-10;
    }
    
    if (h != h) h = 1e-11;
    //std::cerr << "k1 " << x << ' ' << v << ' ' << t << ' ' << h << std::endl;

    

    k1[0] = h*v;
    good_force = good_force && GetForce(x,v,t,efield,bfield,
					index_e, index_b, k1[1]);
    k1[1] = h * k1[1];
    //std::cerr << "k2 " << x << ' ' << v << ' ' << k1[0] << ' ' << k1[1] << ' ' << t << std::endl;

    k2[0] = h*(v+k1[1]*a_rkf45[2][1]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[2][1],
					k2[0]/h,t+h*c_rkf45[2],
					efield,bfield,
					index_e, index_b, k2[1]);
    k2[1] = h * k2[1];
    //std::cerr << "k3" << std::endl;

    k3[0] = h*(v+k1[1]*a_rkf45[3][1]+
	         k2[1]*a_rkf45[3][2]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[3][1]+
					k2[0]*a_rkf45[3][2],
					k3[0]/h,t+h*c_rkf45[3],
					efield,bfield, index_e,
					index_b, k3[1]);
    k3[1] = h * k3[1];
    //std::cerr << "k4" << std::endl;

    k4[0] = h*(v+k1[1]*a_rkf45[4][1]+
	         k2[1]*a_rkf45[4][2]+
	         k3[1]*a_rkf45[4][3]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[4][1]+
					k2[0]*a_rkf45[4][2]+
					k3[0]*a_rkf45[4][3],
					k4[0]/h,t+h*c_rkf45[4],
					efield,bfield, index_e,
					index_b, k4[1]);
    k4[1] = h * k4[1];
    //std::cerr << "k5" << std::endl;

    k5[0] = h*(v+k1[1]*a_rkf45[5][1]+
	         k2[1]*a_rkf45[5][2]+
	         k3[1]*a_rkf45[5][3]+
	         k4[1]*a_rkf45[5][4]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[5][1]+
					k2[0]*a_rkf45[5][2]+
					k3[0]*a_rkf45[5][3]+
					k4[0]*a_rkf45[5][4],
					k5[0]/h,t+h*c_rkf45[5],
					efield,bfield, index_e,
					index_b, k5[1]);
    k5[1] = h * k5[1];
    //std::cerr << "k6" << std::endl;

    k6[0] = h*(v+k1[1]*a_rkf45[6][1]+
	         k2[1]*a_rkf45[6][2]+
	         k3[1]*a_rkf45[6][3]+
	         k4[1]*a_rkf45[6][4]+
	       k5[1]*a_rkf45[6][5]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[6][1]+
					k2[0]*a_rkf45[6][2]+
					k3[0]*a_rkf45[6][3]+
					k4[0]*a_rkf45[6][4]+
					k5[0]*a_rkf45[6][5],
					k6[0]/h,t+h*c_rkf45[6],
					efield,bfield,
					index_e, index_b, k6[1]);
    k6[1] = h * k6[1];

    v_1 = v+(k1[1]*b_rkf451[1]+
	     k2[1]*b_rkf451[2]+
	     k3[1]*b_rkf451[3]+
	     k4[1]*b_rkf451[4]+
	     k5[1]*b_rkf451[5]);

    v_2 = v+(k1[1]*b_rkf452[1]+
	     k2[1]*b_rkf452[2]+
	     k3[1]*b_rkf452[3]+
	     k4[1]*b_rkf452[4]+
	     k5[1]*b_rkf452[5]+
	     k6[1]*b_rkf452[6]);

    x_1 = x+(k1[0]*b_rkf451[1]+
	     k2[0]*b_rkf451[2]+
	     k3[0]*b_rkf451[3]+
	     k4[0]*b_rkf451[4]+
	     k5[0]*b_rkf451[5]);

    x_2 = x+(k1[0]*b_rkf452[1]+
	     k2[0]*b_rkf452[2]+
	     k3[0]*b_rkf452[3]+
	     k4[0]*b_rkf452[4]+
	     k5[0]*b_rkf452[5]+
	     k6[0]*b_rkf452[6]);

    dx = x_1 - x_2;
    dv = v_1 - v_2;
    
    dxlen = sqrt(dx.dot(dx));
    dvlen = sqrt(dv.dot(dv));
    //std::cerr << x << ' ' << x_1 << std::endl;
    
    // std::cout << rkf45_stepsize << ' ' <<
    //   dvlen << ' ' << dxlen << std::endl;
    //std::cout << "finished all" << std::endl;
    //  if(count > 3) {
    // //   h = h/10.0;
    //    std::cerr << count << ' ' << dxlen << ' ' << h << std::endl;
    //  }
  }while ((dxlen > CTolerance || dxlen < CTolerance/50.0) && count < 4);
  
  // if(count != 1) 
  //   std::cerr << "change stepsize " << time << ' ' << count << ' ' << h 
  // 	      << ' ' << search_radiusb << ' ' << search_radiuse << std::endl;

  stepsize = h;

  position = x_2;
  velocity = v_2;

  time += h;

  if(dxlen*30 < CTolerance && dxlen != 0) {
      // std::cerr << h << ' ' << pow(CTolerance / dxlen, 0.25) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      // 		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.25)
      // 		<< std::endl;
      //h = 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2);
    //h = 0.840896 * h * pow(CTolerance / dxlen, 0.25);
  }
  stepsize = h;
  //std::cout << "finished all" << std::endl;
  return good_force;
}

bool ChargedParticle::Prop_DT_RKF45_GC(const EFieldMap & efield,
				    const EFieldMap & bfield, 
				    const my_kd_tree_t & index_e,
				    const my_kd_tree_t & index_b,
				    double & stepsize)
{
  Vector3 k1[2],k2[2],k3[2],k4[2],k5[2],k6[2];
  Vector3 v = velocity;
  Vector3 x = position;
  Vector3 v_1,v_2;
  Vector3 x_1,x_2;
  double h = stepsize;
  double dvlen, dxlen;
  double t = 0;
  dvlen = dxlen = CTolerance;
  int count = 0;
  
  Vector3 dx = x_1 - x_2;
  Vector3 dv = v_1 - v_2;
  bool good_force = true;

  do {
    count++;
    if(count == 4) {
      h=h;
    }
    else if(dxlen > CTolerance) {
      // std::cerr << h << ' ' << dxlen << ' ' << pow(CTolerance / dxlen, 0.2) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      //  		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2) << " s"
      //  		<< std::endl;
      h = 0.840896 * h * pow(CTolerance / abs(dxlen), 0.2);
      if(count == 4 || h != h) h = 1e-12;
    }
    else if(dxlen < CTolerance/50.0) {
      // std::cerr << h << ' ' << dxlen << ' ' << pow(CTolerance / dxlen, 0.2) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      //  		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2) << " b"
      //  		<< std::endl;
      h = 0.840896 * h * pow(CTolerance / abs(dxlen), 0.25);
      if(count == 4 || h != h || h > 1e-10 || dxlen == 0) h = 1e-10;
    }
    
    if (h != h) h = 1e-11;
    //std::cerr << "k1 " << x << ' ' << v << ' ' << t << ' ' << h << std::endl;

    

    k1[0] = h*v;
    good_force = good_force && GetForce(x,v,t,efield,bfield,
					index_e, index_b, k1[1]);
    k1[1] = h * k1[1];
    //std::cerr << "k2 " << x << ' ' << v << ' ' << k1[0] << ' ' << k1[1] << ' ' << t << std::endl;

    k2[0] = h*(v+k1[1]*a_rkf45[2][1]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[2][1],
					k2[0]/h,t+h*c_rkf45[2],
					efield,bfield,
					index_e, index_b, k2[1]);
    k2[1] = h * k2[1];
    //std::cerr << "k3" << std::endl;

    k3[0] = h*(v+k1[1]*a_rkf45[3][1]+
	         k2[1]*a_rkf45[3][2]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[3][1]+
					k2[0]*a_rkf45[3][2],
					k3[0]/h,t+h*c_rkf45[3],
					efield,bfield, index_e,
					index_b, k3[1]);
    k3[1] = h * k3[1];
    //std::cerr << "k4" << std::endl;

    k4[0] = h*(v+k1[1]*a_rkf45[4][1]+
	         k2[1]*a_rkf45[4][2]+
	         k3[1]*a_rkf45[4][3]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[4][1]+
					k2[0]*a_rkf45[4][2]+
					k3[0]*a_rkf45[4][3],
					k4[0]/h,t+h*c_rkf45[4],
					efield,bfield, index_e,
					index_b, k4[1]);
    k4[1] = h * k4[1];
    //std::cerr << "k5" << std::endl;

    k5[0] = h*(v+k1[1]*a_rkf45[5][1]+
	         k2[1]*a_rkf45[5][2]+
	         k3[1]*a_rkf45[5][3]+
	         k4[1]*a_rkf45[5][4]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[5][1]+
					k2[0]*a_rkf45[5][2]+
					k3[0]*a_rkf45[5][3]+
					k4[0]*a_rkf45[5][4],
					k5[0]/h,t+h*c_rkf45[5],
					efield,bfield, index_e,
					index_b, k5[1]);
    k5[1] = h * k5[1];
    //std::cerr << "k6" << std::endl;

    k6[0] = h*(v+k1[1]*a_rkf45[6][1]+
	         k2[1]*a_rkf45[6][2]+
	         k3[1]*a_rkf45[6][3]+
	         k4[1]*a_rkf45[6][4]+
	       k5[1]*a_rkf45[6][5]);
    good_force = good_force && GetForce(x+k1[0]*a_rkf45[6][1]+
					k2[0]*a_rkf45[6][2]+
					k3[0]*a_rkf45[6][3]+
					k4[0]*a_rkf45[6][4]+
					k5[0]*a_rkf45[6][5],
					k6[0]/h,t+h*c_rkf45[6],
					efield,bfield,
					index_e, index_b, k6[1]);
    k6[1] = h * k6[1];

    v_1 = v+(k1[1]*b_rkf451[1]+
	     k2[1]*b_rkf451[2]+
	     k3[1]*b_rkf451[3]+
	     k4[1]*b_rkf451[4]+
	     k5[1]*b_rkf451[5]);

    v_2 = v+(k1[1]*b_rkf452[1]+
	     k2[1]*b_rkf452[2]+
	     k3[1]*b_rkf452[3]+
	     k4[1]*b_rkf452[4]+
	     k5[1]*b_rkf452[5]+
	     k6[1]*b_rkf452[6]);

    x_1 = x+(k1[0]*b_rkf451[1]+
	     k2[0]*b_rkf451[2]+
	     k3[0]*b_rkf451[3]+
	     k4[0]*b_rkf451[4]+
	     k5[0]*b_rkf451[5]);

    x_2 = x+(k1[0]*b_rkf452[1]+
	     k2[0]*b_rkf452[2]+
	     k3[0]*b_rkf452[3]+
	     k4[0]*b_rkf452[4]+
	     k5[0]*b_rkf452[5]+
	     k6[0]*b_rkf452[6]);

    dx = x_1 - x_2;
    dv = v_1 - v_2;
    
    dxlen = sqrt(dx.dot(dx));
    dvlen = sqrt(dv.dot(dv));
    //std::cerr << x << ' ' << x_1 << std::endl;
    
    // std::cout << rkf45_stepsize << ' ' <<
    //   dvlen << ' ' << dxlen << std::endl;
    //std::cout << "finished all" << std::endl;
    //  if(count > 3) {
    // //   h = h/10.0;
    //    std::cerr << count << ' ' << dxlen << ' ' << h << std::endl;
    //  }
  }while ((dxlen > CTolerance || dxlen < CTolerance/50.0) && count < 4);
  
  // if(count != 1) 
  //   std::cerr << "change stepsize " << time << ' ' << count << ' ' << h 
  // 	      << ' ' << search_radiusb << ' ' << search_radiuse << std::endl;

  stepsize = h;

  position = x_2;
  velocity = v_2;

  time += h;

  if(dxlen*30 < CTolerance && dxlen != 0) {
      // std::cerr << h << ' ' << pow(CTolerance / dxlen, 0.25) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      // 		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.25)
      // 		<< std::endl;
      //h = 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2);
    //h = 0.840896 * h * pow(CTolerance / dxlen, 0.25);
  }
  stepsize = h;
  //std::cout << "finished all" << std::endl;
  return good_force;
}

bool ChargedParticle::Prop_DT_RKKC45(const EFieldMap & efield,
				     const EFieldMap & bfield,
				     const my_kd_tree_t & index_e,
				     const my_kd_tree_t & index_b,
				     double & stepsize)
{
  Vector3 k1[2],k2[2],k3[2],k4[2],k5[2],k6[2];
  Vector3 v = velocity;
  Vector3 x = position;
  Vector3 v_1,v_2;
  Vector3 x_1,x_2;
  double h = stepsize;
  double dvlen, dxlen;
  double t = 0;
  dvlen = dxlen = CTolerance;
  int count = 0;
  
  Vector3 dx = x_1 - x_2;
  Vector3 dv = v_1 - v_2;
  bool good_force = true;

  do {
    count++;
    if(count == 4) {
      h=h;
    }
    else if(dxlen > CTolerance) {
      // std::cerr << h << ' ' << dxlen << ' ' << pow(CTolerance / dxlen, 0.2) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      //  		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2) << " s"
      //  		<< std::endl;
      h = 0.840896 * h * pow(CTolerance / abs(dxlen), 0.2);
      if(count == 4 || h != h || h < 1e-12) h = 1e-12;
    }
    else if(dxlen < CTolerance/50.0) {
      // std::cerr << h << ' ' << dxlen << ' ' << pow(CTolerance / dxlen, 0.2) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      //  		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2) << " b"
      //  		<< std::endl;
      h = 0.840896 * h * pow(CTolerance / abs(dxlen), 0.25);
      if(count == 4 || h != h || h > 1e-10 || dxlen == 0) h = 1e-10;
    }
    
    if (h != h) h = 1e-11;
    //std::cerr << "k1 " << x << ' ' << v << ' ' << t << ' ' << h << std::endl;

    k1[0] = h*v;
    good_force = good_force && GetForce(x,v,t,efield,bfield,
					index_e, index_b, k1[1]);
    k1[1] = h * k1[1];
    //std::cerr << "k2 " << x << ' ' << v << ' ' << k1[0] << ' ' << k1[1] << ' ' << t << std::endl;

    k2[0] = h*(v+k1[1]*a_rkkc45[2][1]);
    good_force = good_force && GetForce(x+k1[0]*a_rkkc45[2][1],
					k2[0]/h,t+h*c_rkkc45[2],
					efield,bfield,
					index_e, index_b, k2[1]);
    k2[1] = h * k2[1];
    //std::cerr << "k3" << std::endl;

    k3[0] = h*(v+k1[1]*a_rkkc45[3][1]+
	         k2[1]*a_rkkc45[3][2]);
    good_force = good_force && GetForce(x+k1[0]*a_rkkc45[3][1]+
					k2[0]*a_rkkc45[3][2],
					k3[0]/h,t+h*c_rkkc45[3],
					efield,bfield, index_e,
					index_b, k3[1]);
    k3[1] = h * k3[1];
    //std::cerr << "k4" << std::endl;

    k4[0] = h*(v+k1[1]*a_rkkc45[4][1]+
	         k2[1]*a_rkkc45[4][2]+
	         k3[1]*a_rkkc45[4][3]);
    good_force = good_force && GetForce(x+k1[0]*a_rkkc45[4][1]+
					k2[0]*a_rkkc45[4][2]+
					k3[0]*a_rkkc45[4][3],
					k4[0]/h,t+h*c_rkkc45[4],
					efield,bfield, index_e,
					index_b, k4[1]);
    k4[1] = h * k4[1];
    //std::cerr << "k5" << std::endl;

    k5[0] = h*(v+k1[1]*a_rkkc45[5][1]+
	         k2[1]*a_rkkc45[5][2]+
	         k3[1]*a_rkkc45[5][3]+
	         k4[1]*a_rkkc45[5][4]);
    good_force = good_force && GetForce(x+k1[0]*a_rkkc45[5][1]+
					k2[0]*a_rkkc45[5][2]+
					k3[0]*a_rkkc45[5][3]+
					k4[0]*a_rkkc45[5][4],
					k5[0]/h,t+h*c_rkkc45[5],
					efield,bfield, index_e,
					index_b, k5[1]);
    k5[1] = h * k5[1];
    //std::cerr << "k6" << std::endl;

    k6[0] = h*(v+k1[1]*a_rkkc45[6][1]+
	         k2[1]*a_rkkc45[6][2]+
	         k3[1]*a_rkkc45[6][3]+
	         k4[1]*a_rkkc45[6][4]+
	       k5[1]*a_rkkc45[6][5]);
    good_force = good_force && GetForce(x+k1[0]*a_rkkc45[6][1]+
					k2[0]*a_rkkc45[6][2]+
					k3[0]*a_rkkc45[6][3]+
					k4[0]*a_rkkc45[6][4]+
					k5[0]*a_rkkc45[6][5],
					k6[0]/h,t+h*c_rkkc45[6],
					efield,bfield,
					index_e, index_b, k6[1]);
    k6[1] = h * k6[1];

    v_1 = v+(k1[1]*b_rkkc451[1]+
	     k2[1]*b_rkkc451[2]+
	     k3[1]*b_rkkc451[3]+
	     k4[1]*b_rkkc451[4]+
	     k5[1]*b_rkkc451[5]);

    v_2 = v+(k1[1]*b_rkkc452[1]+
	     k2[1]*b_rkkc452[2]+
	     k3[1]*b_rkkc452[3]+
	     k4[1]*b_rkkc452[4]+
	     k5[1]*b_rkkc452[5]+
	     k6[1]*b_rkkc452[6]);

    x_1 = x+(k1[0]*b_rkkc451[1]+
	     k2[0]*b_rkkc451[2]+
	     k3[0]*b_rkkc451[3]+
	     k4[0]*b_rkkc451[4]+
	     k5[0]*b_rkkc451[5]);

    x_2 = x+(k1[0]*b_rkkc452[1]+
	     k2[0]*b_rkkc452[2]+
	     k3[0]*b_rkkc452[3]+
	     k4[0]*b_rkkc452[4]+
	     k5[0]*b_rkkc452[5]+
	     k6[0]*b_rkkc452[6]);

    dx = x_1 - x_2;
    dv = v_1 - v_2;
    
    dxlen = sqrt(dx.dot(dx));
    dvlen = sqrt(dv.dot(dv));
    //std::cerr << x << ' ' << x_1 << std::endl;
    
    // std::cout << rkkc45_stepsize << ' ' <<
    //   dvlen << ' ' << dxlen << std::endl;
    //std::cout << "finished all" << std::endl;
    //  if(count > 3) {
    // //   h = h/10.0;
    //    std::cerr << count << ' ' << dxlen << ' ' << h << std::endl;
    //  }
  }while ((dxlen > CTolerance || dxlen < CTolerance/50.0) && count < 4);
  
  if(count != 1) 
    std::cerr << "change stepsize " << time << ' ' << count << ' ' << h 
  	      << ' ' << search_radiusb << ' ' << dxlen << std::endl;

  stepsize = h;

  position = x_2;
  velocity = v_2;

  time += h;

  if(dxlen*30 < CTolerance && dxlen != 0) {
      // std::cerr << h << ' ' << pow(CTolerance / dxlen, 0.25) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      // 		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.25)
      // 		<< std::endl;
      //h = 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2);
    //h = 0.840896 * h * pow(CTolerance / dxlen, 0.25);
  }
  stepsize = h;
  //std::cout << "finished all" << std::endl;
  return good_force;


}

bool ChargedParticle::Prop_DT_VelVer(const EFieldMap & efield,
				     const EFieldMap & bfield,
				     const my_kd_tree_t & index_e,
				     const my_kd_tree_t & index_b,
				     double & stepsize)
{
  // maclachlan or something
  Vector3 v = velocity;
  Vector3 x = position;
  Vector3 force,v1,x2,x1,B,E;
  double h = stepsize;
  bool good_force = true;

  x1 = x + h * v * 0.5;

  bool goodE = efield.Interpolate_Field(x1, index_e, search_radiuse, E);
  bool goodB = bfield.Interpolate_Field(x1, index_b, search_radiusb, B);
  //force = (E + vel.cross(B)) * c2m_r;

  double bmag = B.magnitude();
  long double omega, tanTheta05, sinTheta, cosTheta, S,C = 0.;
  
  omega = bmag * c2m_r;
  tanTheta05 = tan(0.5*omega*h);
  double tansq = tanTheta05 * tanTheta05;
  sinTheta = 2*tanTheta05 / (1+tansq);
  cosTheta = (1-tansq) / (1+tansq);

  Vector3 vpara = v.dot(B) * B / (bmag*bmag);
  Vector3 vperp = v - vpara;

  Vector3 Epara = c2m_r * E.dot(B) * B / (bmag*bmag);
  Vector3 Eperp = E - Epara;

  double Eperpmag = Eperp.magnitude();

  Vector3 vB = vpara + cosTheta * vperp + sinTheta * (vperp.cross(B)/bmag);
  Vector3 vE = - h * Epara + ((1-cosTheta)*(Eperp.cross(B)) - sinTheta * Eperp) / omega;

  v1 = vB + vE;
  
  x2 = x1 + h * v1 * 0.5;

  velocity = v1;
  position = x2;
  
  time += h;

  stepsize = h;
  

  // if(time == 0) {
  //   good_force = good_force && GetForce(x,v,t,efield,bfield,
  // 					index_e, index_b, oldaccel);
  //   Prop_DT_RKF451(efield, bfield, index_e, index_b, stepsize);
  // }

  // Vector3 v = velocity;
  // Vector3 x = position;
  // Vector3 force;
  // double h = stepsize;
  // good_force = good_force && GetForce(x,v,t,efield,bfield,
  // 				      index_e, index_b, force);

  // Vector3 x_pdt = x + h * v + (4 * force - oldaccel) * h * h / 6.0;

  return good_force;
}


bool ChargedParticle::Prop_DT_RKF452(const EFieldMap & efield,
				    const EFieldMap & bfield, 
				    const my_kd_tree_t & index_e,
				    const my_kd_tree_t & index_b,
				    double & stepsize)
{
  double h = stepsize;
  Vector3 v = velocity;
  Vector3 x = position;
  bool good_force = true;
  if (time == 0) {
    Vector3 force;
    bool good_force = good_force && GetForce(x,v,time,efield,bfield,
					index_e, index_b, force);
    
    //std::cerr << force << std::endl;
    velocity = v - force * 0.5 * h;
  }
  Vector3 E, B;

  bool goodE = efield.Interpolate_Field(position, index_e, search_radiuse, E);
  bool goodB = bfield.Interpolate_Field(position, index_b, search_radiusb, B);

  //std::cerr << position << ' ' << E << ' ' << B << ' ' << h << std::endl;

  //E = GetFakeE(x);
  //  B = GetFakeB(x);

  //oldvelocity = velocity;
  //oldB = B;

  oldposition = position;
  
  UpdateVelocityBoris(h, E, B);
  
  position = x + velocity * h;

  time += h;
  
  return good_force;
	    
}

// from http://www.particleincell.com/2011/vxb-rotation/
void ChargedParticle::UpdateVelocityBoris(double dt, Vector3 E, Vector3 B)
{
  Vector3 v_minus;
  Vector3 v_prime;
  Vector3 v_plus;
	
  Vector3 t;
  Vector3 s;
  double t_mag2;
  int dim;
	
  /*t vector*/
  t = c2m_r * B * 0.5 * dt;
	
  /*magnitude of t, squared*/
  t_mag2 = t.dot(t);
	
  /*s vector*/
  s = 2.0 * t / (1 + t_mag2);

  /*v minus*/
  v_minus = velocity + c2m_r * E * 0.5 * dt;
	
  /*v prime*/
  v_prime = v_minus + v_minus.cross(t);
	
  /*v prime*/
  v_plus = v_minus + v_prime.cross(s);
	
  /*v n+1/2*/
  velocity = v_plus + c2m_r * E * 0.5 * dt;
  vold_half = velocity;
}

bool ChargedParticle::Prop_DT_DOPRI(const EFieldMap & efield,
				    const EFieldMap & bfield, 
				    const my_kd_tree_t & index_e,
				    const my_kd_tree_t & index_b,
				    double & stepsize)
{
  /*
  Vector3 k1[2],k2[2],k3[2],k4[2],k5[2],k6[2],k7[2],k8[2],k9[2],k10[2],k11[2],k12[2],k13[2];
  Vector3 v = velocity;
  Vector3 x = position;
  Vector3 v_1,v_2;
  Vector3 x_1,x_2;
  double h = stepsize;
  double dvlen, dxlen;
  double t = 0;
  dvlen = dxlen = CTolerance;
  int count = 0;
  
  Vector3 dx = x_1 - x_2;
  Vector3 dv = v_1 - v_2;

  do {
    count++;
    if(dxlen > CTolerance || dvlen > CTolerance) {
      //h = 0.840896 * sqrt(sqrt(CTolerance * h / (dvlen)));
      h = 0.840896 * h * pow(CTolerance / dxlen, 0.2);
    }

    //h = rkf45_stepsize;
    k1[0] = h*v;
    k1[1] = h*GetForce(x,v,t,efield,bfield, index_e, index_b);
    //std::cout << "finished k1[0] k1[1]" << std::endl;

    k2[0] = h*(v+k1[1]*a_dp[2][1]);
    k2[1] = h*GetForce(x+k1[0]*a_dp[2][1],k2[0]/h,t+h*c_dp[2],efield,bfield, index_e, index_b);

    k3[0] = h*(v+k1[1]*a_dp[3][1]+
	         k2[1]*a_dp[3][2]);
    k3[1] = h*GetForce(x+k1[0]*a_dp[3][1]+
		         k2[0]*a_dp[3][2],
    		       k3[0]/h,t+h*c_dp[3],efield,bfield, index_e, index_b);

    k4[0] = h*(v+k1[1]*a_dp[4][1]+
	         k2[1]*a_dp[4][2]+
	         k3[1]*a_dp[4][3]);
    k4[1] = h*GetForce(x+k1[0]*a_dp[4][1]+
		         k2[0]*a_dp[4][2]+
		         k3[0]*a_dp[4][3],
    		       k4[0]/h,t+h*c_dp[4],efield,bfield, index_e, index_b);

    k5[0] = h*(v+k1[1]*a_dp[5][1]+
	         k2[1]*a_dp[5][2]+
	         k3[1]*a_dp[5][3]+
	         k4[1]*a_dp[5][4]);
    k5[1] = h*GetForce(x+k1[0]*a_dp[5][1]+
		         k2[0]*a_dp[5][2]+
		         k3[0]*a_dp[5][3]+
		         k4[0]*a_dp[5][4],
		       k5[0]/h,t+h*c_dp[5],efield,bfield, index_e, index_b);

    k6[0] = h*(v+k1[1]*a_dp[6][1]+
	         k2[1]*a_dp[6][2]+
	         k3[1]*a_dp[6][3]+
	         k4[1]*a_dp[6][4]+
	         k5[1]*a_dp[6][5]);
    k6[1] = h*GetForce(x+k1[0]*a_dp[6][1]+
		         k2[0]*a_dp[6][2]+
		         k3[0]*a_dp[6][3]+
		         k4[0]*a_dp[6][4]+
		         k5[0]*a_dp[6][5],
    		         k6[0]/h,t+h*c_dp[6],efield,bfield, index_e, index_b);

    k7[0] = h*(v+k1[1]*a_dp[7][1]+
	         k2[1]*a_dp[7][2]+
	         k3[1]*a_dp[7][3]+
	         k4[1]*a_dp[7][4]+
	         k5[1]*a_dp[7][5]+
	         k6[1]*a_dp[7][6]);
    k7[1] = h*GetForce(x+k1[0]*a_dp[7][1]+
		         k2[0]*a_dp[7][2]+
		         k3[0]*a_dp[7][3]+
		         k4[0]*a_dp[7][4]+
		         k5[0]*a_dp[7][5]+
		         k6[0]*a_dp[7][6],
    		       k7[0]/h,t+h*c_dp[7],efield,bfield, index_e, index_b);

    k8[0] = h*(v+k1[1]*a_dp[8][1]+
	         k2[1]*a_dp[8][2]+
	         k3[1]*a_dp[8][3]+
	         k4[1]*a_dp[8][4]+
	         k5[1]*a_dp[8][5]+
	         k6[1]*a_dp[8][6]+
	         k7[1]*a_dp[8][7]);
    k8[1] = h*GetForce(x+k1[0]*a_dp[8][1]+
		         k2[0]*a_dp[8][2]+
		         k3[0]*a_dp[8][3]+
		         k4[0]*a_dp[8][4]+
		         k5[0]*a_dp[8][5]+
		         k6[0]*a_dp[8][6]+
		         k7[0]*a_dp[8][7],
    		       k8[0]/h,t+h*c_dp[8],efield,bfield, index_e, index_b);

    k9[0] = h*(v+k1[1]*a_dp[9][1]+
	         k2[1]*a_dp[9][2]+
	         k3[1]*a_dp[9][3]+
	         k4[1]*a_dp[9][4]+
	         k5[1]*a_dp[9][5]+
	         k6[1]*a_dp[9][6]+
	         k7[1]*a_dp[9][7]+
	         k8[1]*a_dp[9][8]);
    k9[1] = h*GetForce(x+k1[0]*a_dp[9][1]+
		         k2[0]*a_dp[9][2]+
		         k3[0]*a_dp[9][3]+
		         k4[0]*a_dp[9][4]+
		         k5[0]*a_dp[9][5]+
		         k6[0]*a_dp[9][6]+
		         k7[0]*a_dp[9][7]+
		         k8[0]*a_dp[9][8],
    		       k9[0]/h,t+h*c_dp[9],efield,bfield, index_e, index_b);

    k10[0] = h*(v+k1[1]*a_dp[10][1]+
	          k2[1]*a_dp[10][2]+
	          k3[1]*a_dp[10][3]+
	          k4[1]*a_dp[10][4]+
	          k5[1]*a_dp[10][5]+
	          k6[1]*a_dp[10][6]+
	          k7[1]*a_dp[10][7]+
	          k8[1]*a_dp[10][8]+
	          k9[1]*a_dp[10][9]);
    k10[1] = h*GetForce(x+k1[0]*a_dp[10][1]+
		          k2[0]*a_dp[10][2]+
		          k3[0]*a_dp[10][3]+
		          k4[0]*a_dp[10][4]+
		          k5[0]*a_dp[10][5]+
		          k6[0]*a_dp[10][6]+
		          k7[0]*a_dp[10][7]+
		          k8[0]*a_dp[10][8]+
	                  k9[0]*a_dp[10][9],
    		       k10[0]/h,t+h*c_dp[10],efield,bfield, index_e, index_b);

    k11[0] = h*(v+k1[1]*a_dp[11][1]+
	          k2[1]*a_dp[11][2]+
	          k3[1]*a_dp[11][3]+
	          k4[1]*a_dp[11][4]+
	          k5[1]*a_dp[11][5]+
	          k6[1]*a_dp[11][6]+
	          k7[1]*a_dp[11][7]+
	          k8[1]*a_dp[11][8]+
	          k9[1]*a_dp[11][9]+
	         k10[1]*a_dp[11][10]);
    k11[1] = h*GetForce(x+k1[0]*a_dp[11][1]+
		          k2[0]*a_dp[11][2]+
		          k3[0]*a_dp[11][3]+
		          k4[0]*a_dp[11][4]+
		          k5[0]*a_dp[11][5]+
		          k6[0]*a_dp[11][6]+
		          k7[0]*a_dp[11][7]+
		          k8[0]*a_dp[11][8]+
	                  k9[0]*a_dp[11][9]+
	                 k10[0]*a_dp[11][10],
    		       k11[0]/h,t+h*c_dp[11],efield,bfield, index_e, index_b);

    k12[0] = h*(v+k1[1]*a_dp[12][1]+
	          k2[1]*a_dp[12][2]+
	          k3[1]*a_dp[12][3]+
	          k4[1]*a_dp[12][4]+
	          k5[1]*a_dp[12][5]+
	          k6[1]*a_dp[12][6]+
	          k7[1]*a_dp[12][7]+
	          k8[1]*a_dp[12][8]+
	          k9[1]*a_dp[12][9]+
	         k10[1]*a_dp[12][10]+
	         k11[1]*a_dp[12][11]);
    k12[1] = h*GetForce(x+k1[0]*a_dp[12][1]+
		          k2[0]*a_dp[12][2]+
		          k3[0]*a_dp[12][3]+
		          k4[0]*a_dp[12][4]+
		          k5[0]*a_dp[12][5]+
		          k6[0]*a_dp[12][6]+
		          k7[0]*a_dp[12][7]+
		          k8[0]*a_dp[12][8]+
	                  k9[0]*a_dp[12][9]+
	                 k10[0]*a_dp[12][10]+
	                 k11[0]*a_dp[12][11],
    		       k12[0]/h,t+h*c_dp[12],efield,bfield, index_e, index_b);

    k13[0] = h*(v+k1[1]*a_dp[13][1]+
	          k2[1]*a_dp[13][2]+
	          k3[1]*a_dp[13][3]+
	          k4[1]*a_dp[13][4]+
	          k5[1]*a_dp[13][5]+
	          k6[1]*a_dp[13][6]+
	          k7[1]*a_dp[13][7]+
	          k8[1]*a_dp[13][8]+
	          k9[1]*a_dp[13][9]+
	         k10[1]*a_dp[13][10]+
	         k11[1]*a_dp[13][11]+
	         k12[1]*a_dp[13][12]);
    k13[1] = h*GetForce(x+k1[0]*a_dp[13][1]+
		          k2[0]*a_dp[13][2]+
		          k3[0]*a_dp[13][3]+
		          k4[0]*a_dp[13][4]+
		          k5[0]*a_dp[13][5]+
		          k6[0]*a_dp[13][6]+
		          k7[0]*a_dp[13][7]+
		          k8[0]*a_dp[13][8]+
	                  k9[0]*a_dp[13][9]+
	                 k10[0]*a_dp[13][10]+
	                 k11[0]*a_dp[13][11]+
			 k12[0]*a_dp[13][12],
    		       k13[0]/h,t+h*c_dp[13],efield,bfield, index_e, index_b);

    v_1 = (k1[1] * b_dp2[1]+
	   k6[1] * b_dp2[6]+
	   k7[1] * b_dp2[7]+
	   k8[1] * b_dp2[8]+
	   k9[1] * b_dp2[9]+
	   k10[1] * b_dp2[10]+
	   k11[1] * b_dp2[11]+
	   k12[1] * b_dp2[12]) + v;

    v_2 = (k1[1] * b_dp1[1]+
	   k6[1] * b_dp1[6]+
	   k7[1] * b_dp1[7]+
	   k8[1] * b_dp1[8]+
	   k9[1] * b_dp1[9]+
	   k10[1] * b_dp1[10]+
	   k11[1] * b_dp1[11]+
	   k12[1] * b_dp1[12]+
	   k13[1] * b_dp1[13]) + v;

    x_1 = (k1[0] * b_dp2[0]+
	   k6[0] * b_dp2[6]+
	   k7[0] * b_dp2[7]+
	   k8[0] * b_dp2[8]+
	   k9[0] * b_dp2[9]+
	   k10[0] * b_dp2[10]+
	   k11[0] * b_dp2[11]+
	   k12[0] * b_dp2[12]) + x;

    x_2 = (k1[0] * b_dp1[0]+
	   k6[0] * b_dp1[6]+
	   k7[0] * b_dp1[7]+
	   k8[0] * b_dp1[8]+
	   k9[0] * b_dp1[9]+
	   k10[0] * b_dp1[10]+
	   k11[0] * b_dp1[11]+
	   k12[0] * b_dp1[12]+
	   k13[0] * b_dp1[13]) + x;

    dx = x_1 - x_2;
    dv = v_1 - v_2;
    
    dxlen = sqrt(dx.dot(dx));
    dvlen = sqrt(dv.dot(dv))*h;
    
    // std:: cout << rkf45_stepsize << ' ' <<
    //   dvlen << ' ' << dxlen << std::endl;
    //std::cerr << x_1 << ' ' << x_2  << std::endl;
  }while ((dxlen > CTolerance || dvlen > CTolerance));

  // if(dvlen*10 < CTolerance && dxlen*10 < CTolerance) {
  //   //h = 0.840896 * sqrt(sqrt(CTolerance * h / (dvlen)));
  //   h = 0.840896 * h * pow(CTolerance / dxlen, 0.25);
  // }

  stepsize = h;

  position = x_1;
  velocity = v_1;

  time += h;
  //std::cout << "finished all" << std::endl;
  */
}

void ChargedParticle::Prop_DT_DOPRI2(const EFieldMap & efield,
				    const EFieldMap & bfield, 
				    const my_kd_tree_t & index_e,
				    const my_kd_tree_t & index_b,
				    double & stepsize)
{
  /*
  Vector3 k[13][2];
  Vector3 v = velocity;
  Vector3 x = position;
  Vector3 v_1,v_2;
  Vector3 x_1,x_2;
  double h = stepsize;
  double dvlen, dxlen;
  double t = 0;
  dvlen = dxlen = CTolerance;
  int count = 0;
  
  Vector3 dx = x_1 - x_2;
  Vector3 dv = v_1 - v_2;

  do {
    count++;
    if(dxlen > CTolerance || dvlen > CTolerance) {
      //h = 0.840896 * sqrt(sqrt(CTolerance * h / (dvlen)));
      h = 0.840896 * h * pow(CTolerance / dxlen, 0.2);
    }

    //h = rkf45_stepsize;
    k[0][0] = h*v;
    k[0][1] = h*GetForce(x,v,t,efield,bfield, index_e, index_b);
    //std::cout << "finished k[0][0] k[0][1]" << std::endl;

    for(int i = 1; i <= 5; i++) {
      k[i][0] = v;
      //std::cerr << "v ";
      for(int j = 0; j < i; j++) {
	k[i][0] = k[i][0] + k[j][1]*a_dp[i+1][j+1];
	//std::cerr <<  "+k[" << j << "][1]a[" << i+1 << "][" << j+1 << "]";
      }
      k[i][0] = k[i][0] * h;
      k[i][1] = x;
      //std::cerr << std::endl << "x ";
      
      for(int j = 0; j < i; j++) {
	k[i][1] = k[i][1] + k[j][0]*a_dp[i+1][j+1];
	//std::cerr <<  "+k[" << j << "][0]a[" << i+1 << "][" << j+1 << "]";
      }
      //std::cerr << std::endl;
      k[i][1] = h * GetForce(k[i][1],k[i][0]/h,t+h*c_dp[i+1],
			     efield, bfield, index_e, index_b);
    }
    std::cerr << k[5][0] << std::endl;
    x_1=x; x_2=x;
    v_1=v; v_2=v;
    for(int i = 0; i < 12; i++) {
      x_1 = x_1 + k[i][0]*b_dp2[i+1];
      v_1 = v_1 + k[i][1]*b_dp2[i+1];
    }
    for(int i = 0; i < 13; i++) {
      x_2 = x_2 + k[i][0]*b_dp1[i+1];
      v_2 = v_2 + k[i][1]*b_dp1[i+1];
    }

    dx = x_1 - x_2;
    dv = v_1 - v_2;
    
    dxlen = sqrt(dx.dot(dx));
    dvlen = sqrt(dv.dot(dv))*h;
    
    // std:: cout << rkf45_stepsize << ' ' <<
    //   dvlen << ' ' << dxlen << std::endl;
    //std::cerr << x_1 << ' ' << x_2  << std::endl;
  }while ((dxlen > CTolerance));

  // if(dvlen*10 < CTolerance && dxlen*10 < CTolerance) {
  //   //h = 0.840896 * sqrt(sqrt(CTolerance * h / (dvlen)));
  //   h = 0.840896 * h * pow(CTolerance / dxlen, 0.25);
  // }
  
  if(count != 1) 
    std::cerr << "change stepsize " << time << ' ' << count << ' ' << h << std::endl;

  stepsize = h;

  position = x_1;
  velocity = v_1;

  time += h;

  if(dxlen*10 < CTolerance && dxlen != 0) {
      // std::cerr << h << ' ' << pow(CTolerance / dxlen, 0.25) << ' '
      // 		<< (CTolerance / dxlen) << ' ' 
      // 		<< 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.25)
      // 		<< std::endl;
      //h = 0.840896 * h * pow((CTolerance / abs(dxlen)), 0.2);
    h = 0.840896 * h * pow(CTolerance / dxlen, 0.25);
  }
  stepsize = h;
  */
}

Vector3 ChargedParticle::GetForce(Vector3 point, Vector3 vel, double t,
				  const EFieldMap & efield,
				  const MagFieldMap & bfield, const my_kd_tree_t & index)
{
  Vector3 E(0,0,0);// = efield.Interpolate_Field(point, index, search_radiuse);
  //Vector3 B = bfield.FieldEval(point);
  //std::cout << result << " " << E << ' ' << B << std::endl;
  //Vector3 E(0,0,0);
  Vector3 B(5,0,0);
  return (E + vel.cross(B)) * charge / mass;
  //return (vel.cross(B)) * charge / mass;
}

bool ChargedParticle::GetForce(Vector3 point, Vector3 vel, double t,
				  const EFieldMap & efield,
				  const EFieldMap & bfield,
				  const my_kd_tree_t & index_e,
				  const my_kd_tree_t & index_b,
				  Vector3 & force)
{
  bool nonrel = true;
  Vector3 E, B;
  if(!nonrel && false) {
    double invcsq = 1.0/(299792458.0*299792458.0);
    double invgamsqrt = sqrt(1 - vel.dot(vel) * invcsq);
    bool goodE = efield.Interpolate_Field(point, index_e, search_radiuse, E);
    bool goodB = bfield.Interpolate_Field(point, index_b, search_radiusb, B);
    force = (E + vel.cross(B) - invcsq * vel * (vel.dot(E))) * c2m_r;
    return (goodE && goodB);
  }
  else {
    bool goodE = efield.Interpolate_Field(point, index_e, search_radiuse, E);
    bool goodB = bfield.Interpolate_Field(point, index_b, search_radiusb, B);

    double bmag = B.magnitude();
    
    Vector3 vpara = vel.projection(B);
    Vector3 vperp = vel - vpara;

    Vector3 B_hat = B / B.magnitude();

    double vperpmag = vperp.magnitude();

    Vector3 bforce = vel.cross(B);
    

    std::cerr << point << ' ' << E << ' ' << B << std::endl;
    force = (E + bforce) * c2m_r;
    //Vector3 temp = (bforce.projection(B));
    //std::cerr << temp << std::endl;
    
    return (goodE && goodB);
    //force = (E);

    Vector3 vparahat = vpara / vpara.magnitude();
    Vector3 vperphat = vperp / vperp.magnitude();
    
    Vector3 tforce, tB;
    double radius = vperpmag / (c2m_r * bmag);

    Vector3 rot_axis = point + 
      radius * (vel.cross(B)) / (vel.magnitude() * B.magnitude());

    double s_theta = 1;
    double c_theta = 0;
    c_theta=sqrt2o2;
    s_theta=sqrt2o2;
    double arr[] = {c_theta + B_hat.xc() * B_hat.xc() * (1-c_theta),
	       B_hat.xc()*B_hat.yc()*(1-c_theta)-B_hat.zc()*s_theta,
	       B_hat.xc()*B_hat.zc()*(1-c_theta)+B_hat.yc()*s_theta,
	       B_hat.xc()*B_hat.yc()*(1-c_theta)+B_hat.zc()*s_theta,
	       c_theta + B_hat.yc() * B_hat.yc() *(1-c_theta),
	       B_hat.zc()*B_hat.yc()*(1-c_theta)-B_hat.xc()*s_theta,
	       B_hat.xc()*B_hat.zc()*(1-c_theta)-B_hat.yc()*s_theta,
	       B_hat.zc()*B_hat.yc()*(1-c_theta)+B_hat.xc()*s_theta,
	       c_theta + B_hat.zc() * B_hat.zc() *(1-c_theta)};
    std::vector<double> elems (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    Matrix mat(elems);
    Vector3 tnewpos = point;
    Vector3 tnewvel = vel;
    tforce = (vel.cross(B));
    //std::cerr << t << ' ' << tnewpos << ' ';
    for (int i = 1; i < 8; i++) {
      tnewpos = mat.Multiply(tnewpos);
      tnewvel = mat.Multiply(tnewvel);
      //tnewpos = mat.Multiply(tnewpos);
      tforce = tforce + tnewvel.cross(B);
      //std::cerr << tnewvel << ' ';
      //std::cerr << tnewpos << ' ';
    }
    //std::cerr << std::endl;
    Vector3 tv = vel.cross(B);

    Vector3 t1 = 0.125 * tforce.projection(B_hat) * c2m_r;
    Vector3 t0 = 0.125 * tforce * c2m_r;
    Vector3 t2 = tv * c2m_r;
    Vector3 t3 = c2m_r*(E + 0.125 * tforce.projection(B_hat)
		  + (vel.cross(B)) - (vel.cross(B)).projection(B_hat));
    //std::cerr << t0 << ' ' << t1 << ' ' << t2 << ' ' << t3 << ' ' << std::endl;
    
    //force = (E + 0.125 * tforce.projection(B_hat)
    //	     + (vel.cross(B)) - (vel.cross(B)).projection(B_hat));
    force = (E + 0.125 * tforce.projection(B_hat));
    //force = tforce * 0.125 + E;
    force = force * c2m_r;
    //std::cerr << t << ' ' << force << std::endl;
    return (goodE && goodB);
  }
}

bool ChargedParticle::GetForce_E(Vector3 point, Vector3 vel, double t,
				    const EFieldMap & efield,
				    const my_kd_tree_t & index_e,
				    Vector3 & field)
{
  bool goodE = true;//efield.Interpolate_Field(point, index_e, search_radiuse, field);
  field = GetFakeE(point);
  field = field * c2m_r;
  return goodE;
}

bool ChargedParticle::GetForce_B(Vector3 point, Vector3 vel, double t,
				    const EFieldMap & bfield,
				    const my_kd_tree_t & index_b,
				    Vector3 & field)
{
  bool goodB = true;//bfield.Interpolate_Field(point, index_b, search_radiusb, field);
  field = GetFakeB(point);
  field = (vel.cross(field)) * c2m_r;
  return goodB;
}

bool ChargedParticle::Get_E(Vector3 point, Vector3 vel, double t,
			    const EFieldMap & efield,
			    const my_kd_tree_t & index_e,
			    Vector3 & field)
{
  bool goodE = true;//efield.Interpolate_Field(point, index_e, search_radiuse, field);
  field = GetFakeE(point);

  return goodE;
}

bool ChargedParticle::Get_B(Vector3 point, Vector3 vel, double t,
			    const EFieldMap & bfield,
			    const my_kd_tree_t & index_b,
			    Vector3 & field)
{
  bool goodB = true;//bfield.Interpolate_Field(point, index_b, search_radiusb, field);
  //field = oldB;
  field = GetFakeB(point);
  return goodB;
}

Vector3 ChargedParticle::GetFakeB(Vector3 point)
{

  double phi = atan(point.yc()/point.xc());

  if(point.xc() < 0) {
    phi += pi;
  }

  Vector3 field(-5.0 * 0.3 * sin(phi), 5.0 * 0.3 * cos(phi), 0);


  // Vector3 field = Vector3(5.0, 0, 0);
  // double x = point.xc();
  // if(x < 0.15) {
  //   field = Vector3(-intpow((0.15-x),2) + 5.0, intpow(x-0.15,2), 0);
  //   //field = Vector3(-(x-0.15)*10.0 + 5.0, 0, 0);
  // }
  // else if(x > 0.3) {
  //   field = Vector3(-intpow((x-0.3),2) + 5.0, intpow(x-0.3,2), 0);
  //   //field = Vector3(-(x-0.3)*10.0 + 5.0, 0, 0);
  // }
  // //std::cerr << x << ' ' << field.xc() << std::endl;
  return field;
}
Vector3 ChargedParticle::GetFakeE(Vector3 point)
{
  //return Vector3(0,0,500);
  return NullVector;
}

void ChargedParticle::SetOldPos()
{
  oldposition = position;
}

Vector3 ChargedParticle::GetPos()
{
  return position;
}

Vector3 ChargedParticle::GetOldPos()
{
  return oldposition;
}

Vector3 ChargedParticle::GetOldVel()
{
  return oldvelocity;
}

Vector3 ChargedParticle::GetOldB()
{
  return oldB;
}

Vector3 ChargedParticle::GetVel()
{
  return velocity;
}

double ChargedParticle::GetPE()
{
  return pot_E;
}

double ChargedParticle::GetKE()
{
  return kin_E;
}
double ChargedParticle::GetCharge()
{
  return charge;
}
double ChargedParticle::GetTime()
{
  return time;
}
double ChargedParticle::GetSRB()
{
  return search_radiusb;
}
double ChargedParticle::GetSRE()
{
  return search_radiuse;
}

// ChargedParticle& ChargedParticle::operator=(const ChargedParticle& rhs)
// {
//   if (this == &rhs) 
//     return *this;

//   double tx, ty, tz;
//   rhs.return_coords(tx, ty, tz);
//   set_coords(tx, ty, tz);
  
//   return *this;
// }


// Overloaded friends!!

std::ostream& operator<< (std::ostream &out, ChargedParticle &tChaPar)
{
    // Since operator<< is a friend of the Point class, we can access
    // Point's members directly.
  double time = tChaPar.GetTime();
  Vector3 tpos = tChaPar.GetPos();
  Vector3 tvel = tChaPar.GetVel();
  Vector3 tovel = tChaPar.GetOldVel();
  Vector3 topos = tChaPar.GetOldPos();
  out << time << " " << tpos << " " << tvel << " " << topos;
  return out;
}

// std::istream& operator>> (std::istream &in, ChargedParticle &tVect)
// {
//     // Since operator<< is a friend of the Point class, we can access
//     // Point's members directly.
//   double tx, ty, tz;
//   in >> tx >> ty >> tz;
//   tVect.set_coords(tx, ty, tz);
//   return in;
// }

// ChargedParticle operator+(const ChargedParticle &v1, const ChargedParticle &v2)
// {
//   return ChargedParticle(v1.xc()+v2.xc(), v1.yc()+v2.yc(), v1.zc()+v2.zc());
// }

// ChargedParticle operator+(const ChargedParticle &v1, const double & scalar)
// {
//   return ChargedParticle(v1.xc()+scalar, v1.yc()+scalar, v1.zc()+scalar);
// }

// ChargedParticle operator-(const ChargedParticle &v1, const ChargedParticle &v2)
// {
//   return ChargedParticle(v1.xc()-v2.xc(), v1.yc()-v2.yc(), v1.zc()-v2.zc());
// }

// ChargedParticle operator-(const ChargedParticle &v1, const double & scalar)
// {
//   return ChargedParticle(v1.xc()-scalar, v1.yc()-scalar, v1.zc()-scalar);
// }

// ChargedParticle operator*(const ChargedParticle &v1, const double &scalar)
// {
//   return ChargedParticle(v1.xc()*scalar, v1.yc()*scalar, v1.zc()*scalar);
// }

