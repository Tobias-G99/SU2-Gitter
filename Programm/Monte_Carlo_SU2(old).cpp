/*
WICHTIG:
Alle benutzten Matrizen benutzen den Stil
array[0]: links oben, array[1]: rechts oben, -conj(array[1]): links unten, conj(array[0]): rechts unten

Hinweis:
Die Kompilation sollte mit dem Zusatzbefehl -std=c++14 durchgeführt werden
*/


//BEGINN DER BIBLIOTHEKEN

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex>
#include <random>
#include <ctime>

using namespace std;

//ENDE DER BIBLIOTHEKEN



//BEGINN DER GLOBALEN VARIABLEN

//BEGINN VERÄNDERBARER VARIABLEN

//Anzahl der Gitterpunkte in Zeitrichtung
int N_t;
//Anzahl der Gitterpunkte in Raumrichtungen
int N_s;
//Inverse Eichkopplungskonstante
double beta;
//Maximalen Drehwinkel bei der Erzeugung der SU2-Matrix nahe der Einheitsmatrix. Größerer Wert liefert kleinere Akzeptanz.
//Es sollte eine Akzeptanz zwischen 40%-60% angestrebt werden
double delta;
//Anzahl an ermittelten Werten der Observable
int N_cf;
//Anzahl an übersprungenen berechneten Observablen beim Start, (N_ini*N_cor) Sweeps werden übersprungen
int N_ini;
//Aktualisierungen des Links innerhalb eines Metropolis Sweeps
int N_hit;
//Anzahl an sweeps, bis eine Observable gespeichert wird
int N_cor;
//bin size
int b_s;
//Anzahl an erstellten Bootstraps zur Bestimmung des statistischen Fehlers
int N_bt;
//Seed für den Mersenne Twister Pseudo Number Generator
int seed=5555;

//ENDE VERÄNDERBARE VARIABLEN

//Deklaration von Pi
const double pi = 3.141592653589793238462643383279502884;

//Gibt die Genauigkeit bei der Überprüfung, ob es sich um eine SU(2) Matrix handelt, an
double su2_acc = 1e-12;

//ENDE DER GLOBALEN VARIABLEN



//BEGINN DER DEFINITION VON STRUCTS

//Ein Punkt des Gitters; Enthält Informationen über seine Koordinaten und die 4 Linkvariablen in Vorwärtsrichtung
typedef struct{
	int coord[4];				//Koordinaten des Punktes im Stil t,x,y,z
	complex<double> link[4][2];	//Linkvariablen in Vorwärtsrichtung im Stil t,x,y,z ; Eintrag 0, Eintrag 1
} point;

//Das Gitter; Enthält Informationen über die Länge und die Punkte des Gitters
typedef struct{
	int len[4];					//Länge des Gitters im Stil t,x,y,z
	point***** pts;				//Pointer auf Punkte des Gitters
} lattice;

//ENDE DER DEFINITION VON STRUCTS



//BEGINN DER FUNKTIONEN

/*
Erzeugt eine Liste an Pseudozufallszahlen über den Mersenne Twister Pseudo Number Generator;
Diese sind mit dem gleichen seed reproduzierbar
*/
mt19937 gen(seed);
uniform_real_distribution<double> unif(0.,1.);
/*
Erzeugt eine Zufallszahl im Intervall [a,b)
*/
double uniform(double a, double b){
	return a+unif(gen)*(b-a);
}

/*
Gibt die Determinante der 2x2 Matrix zurück
*/
complex<double> det_matrix(complex<double>* matrix){
	return matrix[0]*conj(matrix[0])+matrix[1]*conj(matrix[1]);
}

/*
Multipliziert zwei 2x2 Matrizen, wobei der dritte übergebene Pointer den Speicherort angibt
*/
complex<double>* matrix_mul(complex<double>* matrix1, complex<double>* matrix2, complex<double>* dest){
	complex<double> a=matrix1[0]*matrix2[0]-matrix1[1]*conj(matrix2[1]);		//u=ac-bd*
	complex<double> b=matrix1[0]*matrix2[1]+matrix1[1]*conj(matrix2[0]); 		//v=ad+bc*
	dest[0]=a; dest[1]=b;
	return dest;
}

/*
Multipliziert zwei 2x2 Matrizen, wobei jedoch die zweite Matrix adjungiert wird; der dritte übergebene Pointer ist der Speicherort
*/
complex<double>* matrix_mul_adj(complex<double>* matrix1, complex<double>* matrix2, complex<double>* dest){
	complex<double> a=matrix1[0]*conj(matrix2[0])+matrix1[1]*conj(matrix2[1]);  //u=ac*+bd*
	complex<double> b=-matrix1[0]*matrix2[1]+matrix1[1]*matrix2[0]; 			//v=-ad+bc
	dest[0]=a; dest[1]=b;
	return dest;
}

/*
Erstellt eine 2x2 Einheitsmatrix
*/
complex<double>* create_identity(){
	complex<double>* matrix = (complex<double>*) malloc(2*sizeof(complex<double>));
	if(matrix==NULL) exit(100);
	matrix[0]=1; matrix[1]=0;
	return matrix;
}

/*
Reskaliert die Matrix zu einer SU(2) Matrix
*/
complex<double>* rescale_su2(complex<double>* matrix){
	double rescale_fac=sqrt(real(matrix[0]*conj(matrix[0]))+real(matrix[1]*conj(matrix[1])));
	matrix[0]=matrix[0]/rescale_fac;
	matrix[1]=matrix[1]/rescale_fac;
	return matrix;
}

/*
Erzeugt eine Drehmatrix in der Nähe der Identität
*/
complex<double>* set_close_to_identity(complex<double>* matrix){
	double alpha=uniform(0.,delta);
	double u=uniform(-1.,1.);
	double theta=uniform(0., 2.*pi);
	matrix[0]=cos(alpha)+1i*sin(alpha)*u;											//cos(a)+isin(a)n_3
	matrix[1]=sin(alpha)*(sqrt(1-u*u)*sin(theta)+1i*sqrt(1-u*u)*cos(theta));		//sin(a)(n_2+in_1)
	rescale_su2(matrix);
	return matrix;
}

/*
Erzeugt eine Drehmatrix mit gleichverteilter Wahrscheinlichkeit
*/
complex<double>* set_to_rotation_matrix(complex<double>* matrix){
	double alpha=uniform(0.,2.*pi);
	double u=uniform(-1.,1.);
	double theta=uniform(0., 2.*pi);
	matrix[0]=cos(alpha)+1i*sin(alpha)*u;											//cos(a)+isin(a)n_3
	matrix[1]=sin(alpha)*(sqrt(1-u*u)*sin(theta)+1i*sqrt(1-u*u)*cos(theta));		//sin(a)(n_2+in_1)
	rescale_su2(matrix);
	return matrix;
}

/*
Prüft, ob es sich um eine SU(2) Matrix handelt
return 0: keine SU(2) Matrix, return 1: SU(2) Matrix
*/
bool check_su2(complex<double>* matrix){
	complex<double> res;
	res = det_matrix(matrix);
	if((abs(real(res)-1)<su2_acc) && (imag(res)<su2_acc)) return true;
	else return false;
}

/*
Gibt die Matrix in der Konsole aus
*/
void print_matrix(complex<double>* matrix){
	cout << matrix[0] << " " << matrix[1] << "\n" << -conj(matrix[1]) << " " << conj(matrix[0]) << "\n";
}

/*
Gibt die Matrix in der Konsole aus und prüft, ob es sich um eine SU(2) Matix handelt
*/
void print_su2(complex<double>* matrix){
	print_matrix(matrix);
	if(check_su2(matrix)==1) cout << "Ist eine SU(2) Matrix" << "\n";
	else cout << "Ist keine SU(2) Matrix" << "\n";
}

/*
Initialisiert die Punkte des Gitters und gibt den Pointer auf dieses als Typ struct lattice zurück
*/
lattice* initialize_lattice(int t_len, int x_len, int y_len, int z_len){
	point***** time_space_pts = (point*****) malloc(t_len*sizeof(point****)); if(time_space_pts==NULL) exit(100);
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<t_len;t++){
		time_space_pts[t] = (point****) malloc(x_len*sizeof(point***)); if(time_space_pts[t]==NULL) exit(100);
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<x_len;x++){
			time_space_pts[t][x] = (point***) malloc(y_len*sizeof(point**)); if(time_space_pts[t][x]==NULL) exit(100);
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<y_len;y++){
				time_space_pts[t][x][y] = (point**) malloc(z_len*sizeof(point*)); if(time_space_pts[t][x][y]==NULL) exit(100);
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<z_len;z++){
					time_space_pts[t][x][y][z] = (point*) malloc(sizeof(point)); if(time_space_pts[t][x][y][z]==NULL) exit(100);
					//Setzt die Koordinaten des Punktes
					time_space_pts[t][x][y][z]->coord[0]=t; time_space_pts[t][x][y][z]->coord[1]=x;
					time_space_pts[t][x][y][z]->coord[2]=y; time_space_pts[t][x][y][z]->coord[3]=z;
					//Setzt die Linkvariablen auf die Einheitsmatrix
					for(int i=0;i<4;i++){
						time_space_pts[t][x][y][z]->link[i][0]=1;
						time_space_pts[t][x][y][z]->link[i][1]=0;
					}
				}
			}
		}
	}
	lattice* lat = (lattice*) malloc(sizeof(lattice));
	lat->len[0]=t_len; lat->len[1]=x_len; lat->len[2]=y_len; lat->len[3]=z_len;
	lat->pts=time_space_pts;
	return lat;
}

/*
Befreit den Speicher von einer Variable des Typs struct lattice
*/
void free_lattice(lattice* lat){
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<lat->len[0];t++){
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<lat->len[1];x++){
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<lat->len[2];y++){
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<lat->len[3];z++){
					free(lat->pts[t][x][y][z]);
				}
				free(lat->pts[t][x][y]);
			}
			free(lat->pts[t][x]);
		}
		free(lat->pts[t]);
	}
	free(lat->pts);
	free(lat);
}

/*
Berechnet Summe der plaquetts in Vorwärtsrichtung eines Punktes
*/
complex<double>* calc_plaquette(lattice* lat, complex<double>* P, int t, int x, int y, int z, int mu, int nu){
	int c[4]={t,x,y,z};		//Koordinate des Punkts
	int tc[4]={t,x,y,z};	//Temporäre Koordinate
	int cp1[4];				//Nächster Punkt in Richtung t,x,y,z
	point***** pts = lat->pts;
	
	//Ermittelt die nächsten Punkte in Vorwärtsrichtung unter Berücksichtigung der periodischen Randbedingungen
	for(int i=0;i<4;i++){
		if(c[i]==lat->len[i]-1){
			cp1[i]=0;
		}
		else{
			cp1[i]=c[i]+1;
		}
	}
	
	//Plaquette in Richtung mu, nu
	P[0]=1;P[1]=0;														//mu, nu in den Koordinaten bezeichnen die Vektoren
	matrix_mul(P,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],P);			//U_mu(x)
	tc[mu]=cp1[mu];														//x+mu
	matrix_mul(P,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],P);			//U_nu(x+mu)
	tc[mu]=c[mu]; tc[nu]=cp1[nu];										//x+nu
	matrix_mul_adj(P,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],P);		//U^+_mu(x+nu)
	tc[nu]=c[nu];														//x
	matrix_mul_adj(P,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],P);		//U^+_nu(x)
	return P;
}

/*
Berechnet den Polyakov-loop Wert für einen Raumpunkt
*/
complex<double>* calc_polyakov_loop(lattice* lat, complex<double>* L,int x, int y, int z){
	L[0]=1; L[1]=0;
	point***** pts = lat->pts;
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<lat->len[0];t++){
		matrix_mul(L,pts[t][x][y][z]->link[0],L);
	}
	return L;
}

/*
Berechnet das Staple eines links
*/
complex<double>* calc_staple(lattice* lat, complex<double>* K, int t, int x, int y, int z, int mu){
	complex<double>* tmp_matrix = create_identity();
	int c[4]={t,x,y,z};			//Koordinate des Punkts
	int tc[4]={t,x,y,z};		//Temporäre Koordinate
	int cm1[4]; int cp1[4];		//Vorheriger und nächster Punkt in Richtung t,x,y,z
	point***** pts = lat->pts;
	
	//Ermittelt nächste Punkte in Vor- und Rückwärtsrichtung unter Berücksichtigung der periodischen Randbedingungen
	for(int i=0;i<4;i++){
		if(c[i]==0){
			cm1[i]=lat->len[i]-1;
			cp1[i]=c[i]+1;
		}
		else if(c[i]==lat->len[i]-1){
			cm1[i]=c[i]-1;
			cp1[i]=0;
		}
		else{
			cm1[i]=c[i]-1;
			cp1[i]=c[i]+1;
		}
	}
	
	//Berechnet die Matrixeinträge des Staple K
	K[0]=0; K[1]=0;
	for(int nu=0; nu<4; nu++){
		if(mu!=nu){
			//Vorwärtsorientierte Schleife
			tmp_matrix[0]=1;tmp_matrix[1]=0;													//mu, nu in den Koordinaten bezeichnen die Vektoren
			tc[mu]=cp1[mu];																		//x+mu
			matrix_mul(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);		//U_nu(x+mu)
			tc[mu]=c[mu]; tc[nu]=cp1[nu];														//x+nu
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],tmp_matrix);	//U^+_mu(x+nu)
			tc[nu]=c[nu];																		//x
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);	//U^+_nu(x)
			K[0]+=tmp_matrix[0];K[1]+=tmp_matrix[1];
			//Rückwärtsorientierte Schleife
			tmp_matrix[0]=1;tmp_matrix[1]=0;
			tc[mu]=cp1[mu]; tc[nu]=cm1[nu];														//x+(mu-nu)
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);	//U^+_nu(x+(mu-nu))
			tc[mu]=c[mu];																		//x-nu
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],tmp_matrix);	//U^+_mu(x-nu)
			matrix_mul(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);		//U_nu(x-nu)
			K[0]+=tmp_matrix[0]; K[1]+=tmp_matrix[1];
			tc[nu]=c[nu];																		//x
		}
	}
	free(tmp_matrix);
	return K;
}

/*
Berechnet die Änderung in der Wirkung über das staple K
*/
double calc_delta_S(complex<double>* U_new, complex<double>* U_old, complex<double>* K){
	complex<double> a=U_new[0]-U_old[0]; complex<double> b=U_new[1]-U_old[1];
	return -beta*real(a*K[0]-b*conj(K[1]));														//-beta/2*Re(Tr((U_new-U_old)*K)); Re(Tr(M))=2*Re(M[0])
}

/*
Berechnet die Wirkung des gesamten Gitters
*/
double calc_action(lattice* lat){
	point***** pts=lat->pts;
	int* len=lat->len;
	double S=0;					//Gesamtwirkung des Gitters
	complex<double>* tmp_matrix=create_identity();
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<len[0];t++){
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<len[1];x++){
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<len[2];y++){
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<len[3];z++){
					//Richtung  mu, nu des Plaquettes
					for(int mu=0;mu<3;mu++){
						for(int nu=0;nu<4;nu++){
							if(mu<nu){
								//plaquette Wert für einen Raumzeitpunkt, aufgespannt durch mu, nu
								calc_plaquette(lat,tmp_matrix,t,x,y,z,mu,nu);
								S+=real(2.*(1.-tmp_matrix[0]));	
							}
						}
					}
				}
			}
		}
	}
	free(tmp_matrix);
	return beta/2.*S;
}

/*
Berechnet die gemittelte Summe der plaquettes des Gitters
*/
double calc_average_plaquette(lattice* lat){
	point***** pts=lat->pts;
	int* len=lat->len;
	double P_av=0;			//average plaquette value
	complex<double>* tmp_matrix=create_identity();
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<len[0];t++){
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<len[1];x++){
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<len[2];y++){
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<len[3];z++){
					//Richtung mu,nu des Plaquettes
					for(int mu=0;mu<3;mu++){
						for(int nu=0;nu<4;nu++){
							if(mu<nu){
								//plaquette Wert für einen Raumzeitpunkt, aufgespannt durch mu, nu
								calc_plaquette(lat,tmp_matrix,t,x,y,z,mu,nu);
								P_av+=real(2.*(tmp_matrix[0]));	
							}
						}
					}
				}
			}
		}
	}
	free(tmp_matrix);
	//Normalisierungsfaktor
	return P_av/(12.*len[0]*len[1]*len[2]*len[3]);
}

/*
Berechnet den gemittelten Polyakov-loop Wert
*/
double calc_average_polyakov_loop(lattice* lat){
	point***** pts=lat->pts;
	int* len=lat->len;
	double L_av=0;			//Polyakov loop
	complex<double>* tmp_matrix=create_identity();
	//Schleife über alle Punkte in x-Richtung
	for(int x=0;x<len[1];x++){
		//Schleife über alle Punkte in y-Richtung
		for(int y=0;y<len[2];y++){
			//Schleife über alle Punkte in z-Richtung
			for(int z=0;z<len[3];z++){
				//Wert des polyakov loops für einen Raumpunkt
				calc_polyakov_loop(lat,tmp_matrix,x,y,z);
				L_av+=real(tmp_matrix[0]);	
			}
		}
	}
	free(tmp_matrix);
	//Normalisierungsfaktor
	return L_av/(len[1]*len[2]*len[3]);
}


/*
Überprüft die Eichinvarianz von Observablen, indem für jeden Raumzeitpunkt eine Drehmatrix V erzeugt wird
und die Links zu U'_mu(x) = V(x) U_mu(x) V^dagger(x+mu) aktualisiert werden
*/
lattice* check_for_gauge_invariance(lattice* lat){
	point***** pts=lat->pts;
	int* len=lat->len;
	complex<double>* tmp_matrix=create_identity();	//temporärer Speicher für eine Matrix
	int c;											//Koordinate für x-mu
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<len[0];t++){
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<len[1];x++){
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<len[2];y++){
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<len[3];z++){
					set_to_rotation_matrix(tmp_matrix);
					//aktualisiert alle Matrizen nach U'_mu(x) = V(x) U_mu(x) und U_mu(x-mu) V^dagger(x) unter Beachtung der Randbedingungen
					for(int mu=0;mu<4;mu++){
						matrix_mul(tmp_matrix,pts[t][x][y][z]->link[mu],pts[t][x][y][z]->link[mu]);
						if(mu==0){
							if(t==0) c=len[0]-1; else c=t-1;
							matrix_mul_adj(pts[c][x][y][z]->link[mu],tmp_matrix,pts[c][x][y][z]->link[mu]);
						}
						else if(mu==1){
							if(x==0) c=len[1]-1; else c=x-1;
							matrix_mul_adj(pts[t][c][y][z]->link[mu],tmp_matrix,pts[t][c][y][z]->link[mu]);
						}
						else if(mu==2){
							if(y==0) c=len[2]-1; else c=y-1;
							matrix_mul_adj(pts[t][x][c][z]->link[mu],tmp_matrix,pts[t][x][c][z]->link[mu]);
						}
						else{
							if(z==0) c=len[3]-1; else c=z-1;
							matrix_mul_adj(pts[t][x][y][c]->link[mu],tmp_matrix,pts[t][x][y][c]->link[mu]);
						}
					}
				}
			}
		}
	}
}

/*
Aktualisiert alle links des Gitter N_hit mal
*/
lattice* Metropolis_sweep(lattice* lat){
	complex<double>* tmp_matrix=create_identity();	//temporärer Speicher für eine Matrix
	complex<double>* K = create_identity();			//Staple zum Punkt in Richtung mu
	int accepted=0; int runs=0;						//Für die Bestimmung der acceptance rate
	double r; double delta_S;
	point***** pts=lat->pts; int* len=lat->len;
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<len[0];t++){
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<len[1];x++){
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<len[2];y++){
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<len[3];z++){
					//Schleife für alle 4 Raumrichtungen
					for(int mu=0;mu<4;mu++){
						calc_staple(lat,K,t,x,y,z,mu);
						//Aktualisiert den aktuellen link N_hit mal
						for(int i=0;i<N_hit;i++){
							//modifiziert den link zu U'_mu(x) = U_mu(x) R
							set_close_to_identity(tmp_matrix);
							matrix_mul(pts[t][x][y][z]->link[mu],tmp_matrix,tmp_matrix);
							//Reskaliert die Matrix zu SU(2) im Fall von Maschinenungenauigkeiten
							rescale_su2(tmp_matrix);
							delta_S=calc_delta_S(tmp_matrix,pts[t][x][y][z]->link[mu],K);
							r=uniform(0.,1.);
							if(r<exp(-delta_S)){
								pts[t][x][y][z]->link[mu][0]=tmp_matrix[0];
								pts[t][x][y][z]->link[mu][1]=tmp_matrix[1];
								accepted++;
							}
							runs++;
						}
					}
				}
			}
		}
	}
	
	//Gibt einen Hinweis aus, falls ein Sweep eine acceptance rate hatte, die nicht zwischen 40% und 60% liegt
	if(((double)accepted/runs<0.4)||((double)accepted/runs>0.6)){
		cout << "Hinweis: Akzeptanz des Sweep betrug " << (double)accepted/runs*100 << "%"<< endl;	
	}
	free(tmp_matrix);
	free(K);
	return lat;
}

/*
Berechnet den Mittelwert einer Observable nach dem Monte-Carlo-Verfahren mit Hilfe von Metropolis Sweeps
*/
double* Monte_Carlo(double (*fnc)(lattice*),double* obs){
	double tmp_bin;
	cout << "Initialisiere Gitter..." << endl;
	lattice* lat=initialize_lattice(N_t,N_s,N_s,N_s);
	cout << "Erfolg" << endl;
	cout << "Thermalisiere Gitter..." << endl;
	
	//Überspringt N_ini*N_cor sweeps
	for(int i=0;i<N_ini;i++){
		//Prozentanzeige
		if(i%((int)ceil(N_ini/10.))==0) cout << (double)i*100./N_ini << "%" << endl;
		for(int j=0;j<N_cor;j++){
			Metropolis_sweep(lat);
		}
	}
	cout << "Erfolg" << endl;
	
	cout << "Berechne Werte der Observablen..." << endl;
	//Berechnung der Observablen nach jeweils N_cf*(N_cor+1) sweeps, wobei diese zu b_s großen Bins gemittelt werden
	for(int i=0;i<N_cf;i++){
		//Prozentanzeige
		if(i%((int)ceil(N_cf/10.))==0) cout << (double)i*100./N_cf << "%" << endl;
		tmp_bin=0;
		//Scleife über die bin size
		for(int j=0;j<b_s;j++){
			for(int k=0;k<N_cor+1;k++){	//Die +1 sind ein Fragment aus einer alten Version,
				Metropolis_sweep(lat);	//mit ihr ändert sich nicht die Aussagekraft der Ergebnisse
			}							//und sie wurde dringelassen, um die Ergebnisse reproduzieren zu können
			tmp_bin+=fnc(lat);
		}
		obs[i]=tmp_bin/b_s;
	}
	cout << "Erfolg" << endl;
	free_lattice(lat);
	return obs;
}

/*
Ermittelt den statistischen Fehler der Observablen über bootstrap copies
*/
double calc_error_with_bootstrap(double* obs){
	int bt;
	double wav_obs[N_bt]={0}, wav_wav_obs=0, d_wav_obs=0;
	//Mittelwert
	for(int i=0;i<N_cf;i++){
		wav_obs[0]+=obs[i];	
	}
	wav_obs[0]=wav_obs[0]/N_cf;
	wav_wav_obs+=wav_obs[0];
	//Bootstrap copies
	for(int i=1;i<N_bt;i++){
		for(int j=0;j<N_cf;j++){
			bt=(int)uniform(0,N_cf);
			wav_obs[i]+=obs[bt];
		}
		wav_obs[i]=wav_obs[i]/N_cf;
		wav_wav_obs+=wav_obs[i];
	}
	wav_wav_obs=wav_wav_obs/N_bt;
	//Fehler des Mittelwerts
	for(int i=0;i<N_bt;i++){
		d_wav_obs+=pow(wav_obs[i]-wav_wav_obs,2);
	}
	return sqrt(d_wav_obs/(N_bt-1));
}

//ENDE DER FUNKTIONEN

//BEGINN DER AUSGABEN

/*
Berechnet den gemittelten Plaquettewert und gibt zusätzlich die Werte der Monte Carlo History aus
*/
void P_av_with_history(){
	//Festlegung der Parameter
	N_t=4;
	N_s=8;
	beta=2.3;
	delta=0.95;
	N_cf=100;
	N_ini=10;
	N_hit=20;
	N_cor=40;
	b_s=10;
	N_bt=100000;
	
	clock_t timer=clock();
	//Deklaration der Variablen zur Berechnung und Ausgabe
	double P_av[N_cf]; double wav_P_av=0; double d_wav_P_av=0;
	
	//Simulation für die Werte des average plaquettes
	Monte_Carlo(calc_average_plaquette, &P_av[0]);
	
	//Mittelwert des average plaquettes
	for(int i=0;i<N_cf;i++){
		wav_P_av+=P_av[i];	
	}
	wav_P_av=wav_P_av/N_cf;
	
	//Fehler des average plaquettes
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	
	//Ausgabe der P_av Werte aus der Simulation
	ofstream file;
	file.open("P_av.txt");
	file << "P_av = " << wav_P_av << " +/- " << d_wav_P_av << endl;
	file << endl;
	file << "number " << "P_av" << endl;
	for(int i=0;i<N_cf;i++){
		file << i+1 << " " << P_av[i] << endl;
	}
	file.close();
	
	//Ausgabe der genutzten Parameter
	file.open("P_av_parameter.txt");
	file << "N_t=" << N_t << endl;
	file << "N_s=" << N_s << endl;
	file << "beta=" << beta << endl;
	file << "delta=" << delta << endl;
	file << "N_cf=" << N_cf << endl;
	file << "N_ini=" << N_ini << endl;
	file << "N_hit=" << N_hit << endl;
	file << "N_cor=" << N_cor << endl;
	file << "b_s=" << b_s << endl;
	file << "N_bt=" << N_bt << endl;
	file << "seed=" << seed << endl;
	file.close();
	cout << "Berechnung des gemittelten Plaquettewert nach " << (clock()-timer)/(double) CLOCKS_PER_SEC <<" Sekunden beendet." << endl;
}

/*
Berechnet den average plaquette value für verschiedene Werte von beta
*/
void P_av_vs_beta(){
	//Festlegung der Parameter
	N_t=4;
	N_s=8;
	beta=0.1;
	double delta_arr[15]={3.14,3.14,3.14,3.14,2.35,1.8,1.5,1.3,1.1,1.0,1.0,0.95,0.9,0.85,0.85};
	N_cf=50;
	N_ini=10;
	N_hit=20;
	N_cor=40;
	b_s=10;
	N_bt=100000;
	
	//Deklaration der Variablen zur Berechnung und Ausgabe
	double P_av[N_cf]; double wav_P_av=0; double d_wav_P_av=0;
	
	//Berechnung des average plaquettes gegen beta
	ofstream file;
	file.open("P_av_vs_beta.txt");
	file << "beta" << " " << "P_av" << " " << "d_P_av" << endl;
	for(int i=0;i<15;i++){
		//Übernahme der zuvor ausprobierten delta Werte für best mögliche acceptance rate
		delta=delta_arr[i];
		clock_t timer=clock();
		wav_P_av=0; d_wav_P_av=0;
		
		//Simulation für die Werte des average plaquettes
		Monte_Carlo(calc_average_plaquette, &P_av[0]);
		
		//Mittelwert des average plaquettes
		for(int i=0;i<N_cf;i++){
			wav_P_av+=P_av[i];	
		}
		wav_P_av=wav_P_av/N_cf;
		
		//Fehler des average plaquettes
		d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
		
		//Ausgabe des average plaquettes gegen beta
		file << beta << " " << wav_P_av << " " << d_wav_P_av << endl;
		cout << "Berechnung des gemittelten Plaquettewert für beta = " << beta << " nach " << (clock()-timer)/(double) CLOCKS_PER_SEC <<" Sekunden beendet." << endl;
		beta+=0.2;
	}
	file.close();
	
	//Ausgabe der genutzten Parameter
	file.open("P_av_vs_beta_parameter.txt");
	file << "N_t=" << N_t << endl;
	file << "N_s=" << N_s << endl;
	file << "N_cf=" << N_cf << endl;
	file << "N_ini=" << N_ini << endl;
	file << "N_hit=" << N_hit << endl;
	file << "N_cor=" << N_cor << endl;
	file << "b_s=" << b_s << endl;
	file << "N_bt=" << N_bt << endl;
	file << "seed=" << seed << endl;
	file.close();
}

/*
Betrachtet die Entwicklung des Fehlers des gemittelten plaquette value
über größere Anzahl an bootstrap copies
*/
void P_av_bootstrap_development(double f_num){
	clock_t timer=clock();
	//Deklaration der Variablen zur Berechnung und Ausgabe
	double P_av[N_cf]; double wav_P_av=0; double d_wav_P_av=0;
	wav_P_av=0; d_wav_P_av=0;
	
	//Simulation für die Werte des average plaquettes
	Monte_Carlo(calc_average_plaquette, &P_av[0]);
	
	//Name der txt Datei
	ofstream file;
	file.open("d_wav_P_av_vs_N_bt_" + to_string(f_num) +".txt");
	
	//Ausgabe für verschiedene Werte N_bt
	file << "N_bt" << " " << "d_wav_P_av" << " " << endl;
	N_bt=10;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	N_bt=100;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	N_bt=1000;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	N_bt=10000;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	N_bt=50000;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	N_bt=100000;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	N_bt=200000;
	d_wav_P_av=calc_error_with_bootstrap(&P_av[0]);
	file << N_bt << " " << d_wav_P_av << endl;
	file.close();
	cout << "Berechnung der Fehler für verschiedene N_bt nach " << (clock()-timer)/(double) CLOCKS_PER_SEC <<" Sekunden beendet." << endl;
}

/*
Betrachtet die Entwicklung des Fehlers des gemittelten plaquette value
über eine größere bin size
*/
void P_av_bootstrap_bin_size(){
	//Festlegung der Parameter
	N_t=4;
	N_s=8;
	beta=2.3;
	delta=0.95;
	N_cf=100;
	N_ini=10;
	N_hit=20;
	N_cor=40;
	b_s=1;
	
	for(int i=0;i<7;i++){
		P_av_bootstrap_development(b_s);
		//Vergrößerung der bin size
		b_s=2*b_s;
	}
	
	//Ausgabe der genutzten Parameter
	ofstream file;
	file.open("d_wav_P_av_vs_N_bt_parameter.txt");
	file << "N_t=" << N_t << endl;
	file << "N_s=" << N_s << endl;
	file << "beta=" << beta << endl;
	file << "delta=" << delta << endl;
	file << "N_cf=" << N_cf << endl;
	file << "N_ini=" << N_ini << endl;
	file << "N_hit=" << N_hit << endl;
	file << "N_cor=" << N_cor << endl;
	file << "seed=" << seed << endl;
	file.close();
}

/*
Betrachtet die Entwicklung des Fehlers des gemittelten plaquette value
über eine größeren N_cor Wert
*/
void P_av_bootstrap_N_cor(){
	//Festlegung der Parameter
	N_t=4;
	N_s=8;
	beta=2.3;
	delta=0.95;
	N_cf=100;
	N_ini=10;
	N_hit=20;
	N_cor=1;
	b_s=10;
	
	for(int i=0;i<20;i++){
		P_av_bootstrap_development(N_cor);
		//Vergrößerung der bin size
		N_cor=2*N_cor;
	}
	
	//Ausgabe der genutzten Parameter
	ofstream file;
	file.open("d_wav_P_av_vs_N_bt_parameter.txt");
	file << "N_t=" << N_t << endl;
	file << "N_s=" << N_s << endl;
	file << "beta=" << beta << endl;
	file << "delta=" << delta << endl;
	file << "N_cf=" << N_cf << endl;
	file << "N_ini=" << N_ini << endl;
	file << "N_hit=" << N_hit << endl;
	file << "b_s=" << b_s << endl;
	file << "seed=" << seed << endl;
	file.close();
}


void L_av_histogram(){
	//Festlegung der Parameter
	N_t=4;
	N_s=8;
	beta=2.3;
	delta=0.95;
	N_cf=2000;
	N_ini=10;
	N_hit=10;
	N_cor=400;
	b_s=1;
	N_bt=100000;
	
	clock_t timer=clock();
	//Deklaration der Variablen zur Berechnung und Ausgabe
	double L_av[N_cf];
	double L_av_abs[N_cf]; double wav_L_av_abs=0; double d_wav_L_av_abs=0;
	double L2_av[N_cf]; double wav_L2_av=0; double d_wav_L2_av=0;
	
	//Simulation für die Werte des polyakov loop
	Monte_Carlo(calc_average_polyakov_loop, &L_av[0]);
	
	//Absolutwert des polyakov loop und Mittelwert
	for(int i=0;i<N_cf;i++){
		L_av_abs[i]=abs(L_av[i]);
		wav_L_av_abs+=L_av_abs[i];	
	}
	wav_L_av_abs=wav_L_av_abs/N_cf;
	
	//Fehler des Absolutwerts
	d_wav_L_av_abs=calc_error_with_bootstrap(&L_av_abs[0]);
	
	//Quadrat des polyakov loop und Mittelwert
	for(int i=0;i<N_cf;i++){
		L2_av[i]=L_av[i]*L_av[i];
		wav_L2_av+=L2_av[i];
	}
	wav_L2_av=wav_L2_av/N_cf;
	
	//Fehler des Quadrats
	d_wav_L2_av=calc_error_with_bootstrap(&L2_av[0]);
	
	//Ausgabe der Werte
	ofstream file;
	file.open("L_av_"+to_string(beta)+".txt");
	//Speichern der Observablen in der txt Datei "L_av_(beta).txt"
	file << "L_av_abs = " << wav_L_av_abs << " +/- " << d_wav_L_av_abs << endl;
	file << "L^2_av = " << wav_L2_av << " +/- " << d_wav_L2_av << endl;
	file << "Chi = " <<  N_t*N_s*N_s*N_s*(wav_L2_av-wav_L_av_abs*wav_L_av_abs) << " +/- " << N_t*N_s*N_s*N_s*sqrt(d_wav_L2_av*d_wav_L2_av+pow(2*wav_L_av_abs*d_wav_L_av_abs,2)) << endl;
	file << endl;
	file << "number " << "L_av" << endl;
	for(int i=0;i<N_cf;i++){
		file << i+1 << " " << L_av[i] << endl;
	}
	file.close();
	
	//Ausgabe der genutzten Parameter
	file.open("L_av_"+to_string(beta)+"_parameter.txt");
	file << "N_t=" << N_t << endl;
	file << "N_s=" << N_s << endl;
	file << "beta=" << beta << endl;
	file << "delta=" << delta << endl;
	file << "N_cf=" << N_cf << endl;
	file << "N_ini=" << N_ini << endl;
	file << "N_hit=" << N_hit << endl;
	file << "N_cor=" << N_cor << endl;
	file << "b_s=" << b_s << endl;
	file << "N_bt=" << N_bt << endl;
	file << "seed=" << seed << endl;
	cout << "Berechnung des gemittelten Polyakov Loops nach " << (clock()-timer)/(double) CLOCKS_PER_SEC <<" Sekunden beendet." << endl;
}

//ENDE AUSGABEN



//BEGINN MAIN

int main(){
	P_av_with_history();
	//P_av_vs_beta();
	//P_av_bootstrap_N_cor();
	//P_av_bootstrap_histogram();
	//L_av_histogram();
	return 0;
}
