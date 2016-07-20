#include <cmath>
#include <fstream>
#include <iostream>
#include <ctime>
#include <random>
#include <cstdlib>
#include <string>
#include "nr.h" // For interpolation and integration
#include "nrutil.h"
using namespace std;

#define PI 3.14159265
#define w1 1./500. // energy of soft photons, fixed
#define r0 2.81794E-13 // classical electron radius
#define N 10000 // number of points for arrays
#define l 80 // number of w2 values used

#define EPS 1.0e-6
#define JMAX 20
#define FUNC(x) ((*func)(x))
#define NR_END 1
#define FREE_ARG char*

float w2 = 1000.; // energy of high energy photons
float E = w2;
float intValues [l][2][N];
float w2Values [l];

float *Vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */ {
    float *v;
    v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
    return v-nl+NR_END;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
    free((FREE_ARG) (v+nl-NR_END));
}

void polint(float xa[], float ya[], int n, float x, float *y, float *dy){
    int i,m,ns=1;
    float den,dif,dift,ho,hp,w;
    float *c,*d;
    
    dif = fabs(x-xa[1]);
    c = Vector(1,n);
    d = Vector(1,n);
    for(i=1;i<=n;i++){
        if((dift=fabs(x-xa[i])) < dif) {
            ns=i;
            dif = dift;
        }
        c[i]=ya[i];
        d[i]=ya[i];
    }
    *y=ya[ns--];
    for (m=1; m<n; m++) {
        for (i=1; i<=n-m; i++) {
            ho=xa[i]-x;
            hp=xa[i+m]-x;
            den = ho-hp;
            if(den == 0.0){
                cout << "Error in routine polint" << endl;
                continue;
            }
            w=c[i+1]-d[i];
            den=w/den;
            d[i]=hp*den;
            c[i]=ho*den;
        }
        *y += (*dy=(2*ns < (n-m)?c[ns+1]:d[ns--]));
    }
    free_vector(d,1,n);
    free_vector(c,1,n);
}

float qgaus(float (*func)(float), float a, float b){
    int j;
    float xr, xm, dx, s;
    static float x[] = {0.0,0.1488743389,0.4333953941,0.6794095682,0.8650633666,0.9739065285};
    static float w[] = {0.0,0.2955242247,0.2692667193,0.2190863625,0.1494513491,0.0666713443};
    
    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;
    for (j=1; j<=5; j++) {
        dx=xr*x[j];
        s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s *= xr;
}

float I(float epsilon){ // compute spectrum of the produced particle
    return PI*r0*r0/(4*w1*w1*pow(w2,3))*(4*E*E/((E-epsilon)*epsilon)*log(4*w1*(E-epsilon)*epsilon/E)-8*w1*E+2*(2*w1*E-1)*E*E/((E-epsilon)*epsilon)-(1-1./(w1*E))*pow(E,4)/((E-epsilon)*(E-epsilon)*epsilon*epsilon));
}

void CM(int index, float Z[], float F[]){ // Cutpoint method
    float eValues[N];
    for (int i=0; i<N; i++) {
       eValues[i] = intValues[index][0][i];
    }
    int m = N-2;
    int p[m+1];
    for (int j=0; j<m; j++) {
        for (int i=1; i<N; i++) {
            if(F[i]>(float(j)/m)) {
                p[j] = eValues[i];
                break;
            }
        }
    }
    p[N-1] = 1.0;
    for (int i=0; i<N; i++) {
        double U = (double) rand()/RAND_MAX;
        int L = int(ceil(m*U));
        int a = int(p[L]);
        while (U > F[a]) {
            a++;
        }
        Z[i] = eValues[a];
    }
}

int main(){
    srand(time(NULL)); // seed of the random number
    ofstream f;
    f.open("histogramdataall.txt"); // open a file
    for(int d=0; d<l; d++){
        w2Values[d] = w2;
        E = w2;
        float eMin = E/2.*(1-sqrt(1-1./(w1*E)));
        float eMax = E/2.*(1+sqrt(1-1./(w1*E)));
        float e [N];
        for(int i=0; i<N; i++){
            e[i] = eMin + i*(eMax-eMin)/(N-1); // set the epsilon values
        }
        float iArray [N];
        for(int i=0; i<N; i++){
            iArray[i] = qgaus(I,eMin,e[i]); // Call the Gaussian quadrature class
            e[i] /= e[N-1]; // normalize epsilon
            intValues[d][0][i] = e[i];
        }
        for(int i=0; i<N; i++) { // normalize
            iArray[i] /= iArray[N-1];
            intValues[d][1][i] = iArray[i];
        }
        w2 += 100;
    }

    // Binary search using our own w2 value
    cout << "Enter a numerical value: " << endl;
    cin >> w2;
    ofstream g;
    g.open("w2inputvalue.txt"); // open a file
    g << w2 << endl;
    g << w1 << endl;
    g << N << endl;
    g.close();
    int a = 0, b = l-1, index = 0, isExact = 0;
    while(b-a > 1){ // Binary search the table created for the input w2 value
        if(w2 < w2Values[a] || w2 > w2Values[b]){
            cout << "Invalid value, try again." << endl;
            index = -1;
            break;
        }
        else if(w2 == w2Values[a]){
            index = a;
            isExact = 1;
            break;
        }
        else if(w2 == w2Values[b]){
            index = b;
            isExact = 1;
            break;
        }
        int mid = int(float(b+a)/2.);
        if(w2 > w2Values[mid]){ // our value was higher than the midpoint
            a = mid;
            index = b;
        }
        else if(w2 < w2Values[mid]){ // value lower than midpoint
            b = mid;
            index = a;
        }
        else{ // exact value
            index = mid;
            isExact = 1;
            break;
        }
    }
    cout << "The index your value was closest to is " << index << endl;
    cout << isExact << endl;
    
    if(isExact == 1){ // easy case where we input an exact value
        float iValues [N];
        for (int i=0; i<N; i++){
            iValues[i] = intValues[index][1][i];
        }
        float Z [N];
        CM(index, Z, iValues);
        for (int i=0; i<N; i++) {
            f << Z[i] << endl;
        }
    }
    
    else if (index != -1){ // we interpolate only if the input value is within range
        float w2Val1, w2Val2, w2Val3, w2Val4;
        float iVal1 [N], iVal2 [N], iVal3[N], iVal4[N];
        if(index > l-4){
            w2Val1 = w2Values[l-4]; // Isolate the correct w2 and gamma values
            w2Val2 = w2Values[l-3];
            w2Val3 = w2Values[l-2];
            w2Val4 = w2Values[l-1];
            for(int i=0; i<N; i++){
                iVal1[i] = intValues[l-4][1][i];
                iVal2[i] = intValues[l-3][1][i];
                iVal3[i] = intValues[l-2][1][i];
                iVal4[i] = intValues[l-1][1][i];
            }
        }
        else if(index < 3){
            w2Val1 = w2Values[0];
            w2Val2 = w2Values[1];
            w2Val3 = w2Values[2];
            w2Val4 = w2Values[3];
            for(int i=0; i<N; i++){
                iVal1[i] = intValues[0][1][i];
                iVal2[i] = intValues[1][1][i];
                iVal3[i] = intValues[2][1][i];
                iVal4[i] = intValues[3][1][i];
            }
        }
        else{
            w2Val1 = w2Values[index-1];
            w2Val2 = w2Values[index];
            w2Val3 = w2Values[index+1];
            w2Val4 = w2Values[index+2];
            for (int i=0; i<N; i++) {
                iVal1 [i] = intValues[index-1][1][i];
                iVal2 [i] = intValues[index][1][i];
                iVal3 [i] = intValues[index+1][1][i];
                iVal4 [i] = intValues[index+2][1][i];
            }
        }
        
        float fArray [N];
        // Interpolation
        fArray[0] = 0.0;
        for (int i=1; i<N; i++) { // loop over all epsilon
            float fw2;
            float dfw2;
            float wVal [4] = {w2Val1, w2Val2, w2Val3, w2Val4};
            float iVal [4] = {iVal1[i], iVal2[i], iVal3[i], iVal4[i]};
            polint(wVal-1, iVal-1, 4, w2, &fw2, &dfw2); // Call the polynomial interpolation function
            fArray[i] = fw2; // store interpolated values
        }
        cout << "success" << endl;
        
        // Cutpoint method using interpolated values
        float Z [N];
        CM(index,Z,fArray);
        for (int i=0; i<N; i++) {
            f << Z[i] << endl;
        }
    }
    f.close();
}
