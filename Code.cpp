#include<bits/stdc++.h>
using namespace std;

const double R = 83.14;
const double Tc_CO2 = 304.1;  // in kelvin
const double Pc_CO2 = 73.83;   // in bar
const double w_CO2 = 0.228;  // acentric factor for CO2
const double Tc_H2O = 647.1;  // in kelvin
const double Pc_H2O = 220.55;  // in bar
const double w_H2O = 0.346;  // acentric factor for water

double rootZ(double P, double T){
    double m = 0.37464 + 1.54226*w_CO2 - 0.26992*pow(w_CO2, 2);
    double a = (0.457236*pow(R, 2)*pow(Tc_CO2, 2)/Pc_CO2)*pow(1+m*(1-sqrt(T/Tc_CO2)), 2);
    double b = 0.077796*R*Tc_CO2/Pc_CO2;
    double A = (a*P)/pow(R*T, 2);
    double B = (b*P)/(R*T);
    double Z = 1;
    double new_Z = 0;
    double func_Z = pow(Z, 3) - (1-B)*pow(Z, 2) + (A-2*B-3*pow(B, 2))*Z - (A*B - pow(B, 2) - pow(B, 3));
    double der_Z = 3*pow(Z, 2) - 2*(1-B)*Z + (A-2*B-3*pow(B, 2));  
    double tolerance = INT_MAX;
    while(tolerance > 0.0001){
        new_Z = Z - func_Z/der_Z;
        tolerance = abs(new_Z - Z);
        Z = new_Z;
        func_Z = pow(Z, 3) - (1-B)*pow(Z, 2) + (A-2*B-3*pow(B, 2))*Z - (A*B - pow(B, 2) - pow(B, 3));
        der_Z = 3*pow(Z, 2) - 2*(1-B)*Z + (A-2*B-3*pow(B, 2));
    }
    return Z;
}

double fugacity_CO2(double Z, double P, double T){
    double m = 0.37464 + 1.54226*w_CO2 - 0.26992*pow(w_CO2, 2);
    double a = (0.457236*pow(R, 2)*pow(Tc_CO2, 2)/Pc_CO2)*pow(1+m*(1-sqrt(T/Tc_CO2)), 2);
    double b = 0.077796*R*Tc_CO2/Pc_CO2;
    double A = (a*P)/pow(R*T, 2);
    double B = (b*P)/(R*T);
    double del1 = 1 + sqrt(2);
    double del2 = 1 - sqrt(2);
    double fuga = (Z-1) - log(Z-B) - (A/(B*(del2 - del1)))*1*log((Z+del2*B)/(Z+del1*B));
    return exp(fuga);
}

double calc_Ps(double T){
    double temp = 1 - T/Tc_H2O;
    double a1 = -7.8595178;
    double a2 = 1.8440825;
    double a3 = -11.786649;
    double a4 = 22.680741;
    double a5 = -15.9618719;
    double a6 = 1.8012250;
    double func = (Tc_H2O/T)*(a1*temp + a2*pow(temp, 1.5) + a3*pow(temp, 3) + a4*pow(temp, 3.5) + a5*pow(temp, 4) + a6*pow(temp, 7.5));
    return exp(func)*Pc_H2O;
}

double densityH2O(double P, double T){
    double O = T - 273.15;
    double V0 = (1 + 18.1597*pow(10, -3)*O)/(0.9998 + 18.2249*pow(10, -3)*O - 7.9222*pow(10, -6)*pow(O, 2) - 55.4485*pow(10, -9)*pow(O, 3) + 149.7562*pow(10, -12)*pow(O, 4) - 393.2952*pow(10, -15)*pow(O, 5));
    double B = 19654.32 + 147.037*O - 2.2155*pow(O, 2) + 1.0478*pow(10, -2)*pow(O, 3) - 2.2789*pow(10, -5)*pow(O, 4);
    double A1 = 3.2891 - 2.391*pow(10, -3)*O + 2.8446*pow(10, -4)*pow(O, 2) - 2.82*pow(10, -6)*pow(O, 3) + 8.477*pow(10, -9)*pow(O, 4);
    double A2 = 6.245*pow(10, -5) - 3.913*pow(10, -6)*O - 3.499*pow(10, -8)*pow(O, 2) + 7.942*pow(10, -10)*pow(O, 3) - 3.299*pow(10, -12)*pow(O, 4);
    double x = V0 - (V0*P/(B + A1*P + A2*pow(P, 2)));
    return 1/x;
}

double fugacity_H2O(double P, double T, double Ps, double rhoWater){
    double fugacity = Ps*exp(18.0152*(P-Ps)/(rhoWater*R*T));
    return fugacity;
}

double delB(double T){
    double t = -5.279063;
    double beta = 6.187967;
    double del_B = t + beta*pow(pow(10, 3)/T, 0.5);
    return del_B;
}

double henryConst(double T, double rhoWater, double fugacityH2O){
    double n = -0.114535;
    double func = (1-n)* log(fugacityH2O) + n*log(R*T*rhoWater/18) + 2*rhoWater*delB(T);
    return exp(func);
}

double gamma(double T, double P, double m){
    double c1 = -0.0652869;
    double c2 = 0.00016790636;
    double c3 = 40.838951;
    double c4 = 0;
    double c5 = 0;
    double c6 = -0.039266518;
    double c7 = 0;
    double c8 = 0.021157167;
    double c9 = 6.5486487*pow(10, -6);
    double c10 = 0;
    double e1 = -0.01144624;
    double e2 = 2.8274958*pow(10, -5);
    double e3 = 0;
    double e4 = 0;
    double e5 = 0;
    double e6 = 0.013980876;
    double e7 = 0;
    double e8 = -0.014349005;
    double e9 = 0;
    double e10 = 0;
    double y = c1 + c2*T + c3/T + c4*P + c5/P + c6*P/T + c7*T/pow(P, 2) + c8*P/(630-T) + c9*T*log(P) + c10*P/pow(T, 2);
    double e = e1 + e2*T + e3/T + e4*P + e5/P + e6*P/T + e7*T/pow(P, 2) + e8*P/(630-T) + e9*T*log(P) + e10*P/pow(T, 2);
    double x = 2*m*y +2*m*m*e;
    return exp(x);
}

double K0H2O(double T){
    double O = T - 273.15;
    double x = -2.209 + 3.097*pow(10, -2)*O - 1.098*pow(10, -4)*pow(O, 2) + 2.048*pow(10, -7)*pow(O, 3);
    return pow(10, x);
}

int main(){
    double P, T, m;
    // int n = 101;
    fstream fin;
    fin.open("input.txt", ios::in);
    fstream fout;
    fout.open("output.txt", ios::out);
    // double P, T, m;
    // double P = 83;
    // double T = 60 + 273.15;
    // double m = 0.0172;
    fout<<setw(3)<<"T"<<setw(13)<<"P"<<setw(14)<<"m"<<setw(15)<<"phii"<<setw(17)<<"hi(bar)"<<setw(12)<<"Ki"<<setw(15)<<"yi"<<setw(15)<<"xi"<<endl;
    // fout<<endl;
    fout<<"---------------------------------------------------------------------------------------------------------------"<<endl;
    // fout<<endl;
    while(!fin.eof()){
        fin>>T>>P>>m;
        double Z = rootZ(P, T);
        double fugacityCO2 = fugacity_CO2(Z, P, T);
        double Ps = calc_Ps(T);
        double rhoWater = densityH2O(P, T);
        double fugacityH2O = fugacity_H2O(P, T, Ps, rhoWater);
        double h = henryConst(T, rhoWater, fugacityH2O);
        double Gamma = gamma(T, P, m);
        double K_CO2 = h*Gamma/(P*fugacityCO2);
        double K_H2O = (K0H2O(T)/(fugacityH2O*P))*exp((P-1)*18.18/(R*T));
        double yH2O = (1-1/K_CO2)/(1/K_H2O - 1/K_CO2);
        double yn = 1/(1+yH2O);
        double x = yn/K_CO2;

        // fout<<"Z = "<<Z<<endl;
        // fout<<"fugacityco2 = "<<fugacityCO2<<endl;
        // fout<<"fugacityh2o = "<<fugacityH2O<<endl;
        // fout<<"Ps = "<<Ps<<endl;
        // fout<<"rhowater = "<<rhoWater<<endl;
        // fout<<"henryconst = "<<h<<endl;
        // fout<<"KCO2 = "<<K_CO2<<endl;
        // fout<<"Kh2o = "<<K_H2O<<endl;
        // fout<<"yh2o = "<<yH2O<<endl;
        // fout<<"yn = "<<yn<<endl;
        // fout<<"x = "<<x*100<<endl;
        // fout<<"gamma = "<<Gamma<<endl;
        // fout<<"K0H2O = "<<K0H2O(T)<<endl;
        fout<<T<<setw(11)<<P<<setw(15)<<m<<setw(15)<<fugacityCO2<<setw(15)<<h<<setw(15)<<K_H2O<<setw(15)<<yH2O<<setw(15)<<x<<endl;

    }

    fin.close();

    

    
    // fkout<<"20"<<setw(20)<<"P"<<setw(20)<<"Og x"<<setw(20)<<"m"<<setw(20)<<"Calculated x"<<endl;
   
    fout.close();
}