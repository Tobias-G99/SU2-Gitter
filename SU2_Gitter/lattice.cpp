/*
WICHTIG:
Alle benutzten Matrizen benutzen den Stil
array[0]: links oben, array[1]: rechts oben, -conj(array[1]): links unten, conj(array[0]): rechts unten

Hinweis:
Die Kompilation sollte mit dem Zusatzbefehl -std=c++14 durchgeführt werden
*/


//BEGINN DER BIBLIOTHEKEN

#include <iostream>
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

//Anzahl der Gitterpunkte in die Zeitrichtung
int N_t=10;
//Anzahl der Gitterpunkte in Raumrichtungen
int N_s=10;
//Maximalen Drehwinkel bei der Erzeugung der SU2-Matrix nahe der Einheitsmatrix. Größerer Wert liefert kleinere Akzeptanz. Es sollte eine Akzeptanz zwischen 40%-60% angestrebt werden
double delta=0.25;
//Inverse Eichkopplungskonstante
double beta=2.3;
//Aktualisierungen des Links innerhalb eines Metropolis Sweeps
int N_hit=10;
//Anzahl an ermittelten Werten der Observable
int N_cf=50;
//Anzahl der übersprungenen Sweeps bevor ein Wert zur Observable gespeichert wird
int N_cor=0;
//Anzahl an übersprungenen berechneten Observablen beim Start (N_ini*(N_cf+1) Sweeps werden übersprungen)
int N_ini=0;

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
Erzeugt eine Pseudozufallszahl im Intervall [a,b) ; aktuell bei jedem Durchlauf des Programms gleich
*/
random_device rd;
mt19937 gen(rd());
uniform_real_distribution<double> unif(0.,1.);
double uniform(double a, double b){
	return a+unif(rd)*(b-a);
}

/*
Gibt die Determinante der 2x2 Matrix zurück
*/
complex<double> det_matrix(complex<double>* matrix){
	return matrix[0]*conj(matrix[0])+matrix[1]*conj(matrix[1]);
}

/*
Multipliziert zwei 2x2 Matrizen, wobei der dritte übergebene Pointer den Zielort angibt
*/
complex<double>* matrix_mul(complex<double>* matrix1, complex<double>* matrix2, complex<double>* dest){
	complex<double> a=matrix1[0]*matrix2[0]-matrix1[1]*conj(matrix2[1]); //u=ac-bd*
	complex<double> b=matrix1[0]*matrix2[1]+matrix1[1]*conj(matrix2[0]); //v=ad+bc*
	dest[0]=a; dest[1]=b;
	return dest;
}

/*
Multipliziert zwei 2x2 Matrizen, wobei jedoch die zweite Matrix adjungiert wird
*/
complex<double>* matrix_mul_adj(complex<double>* matrix1, complex<double>* matrix2, complex<double>* dest){
	complex<double> a=matrix1[0]*conj(matrix2[0])+matrix1[1]*conj(matrix2[1]); //u=ac*+bd*
	complex<double> b=-matrix1[0]*matrix2[1]+matrix1[1]*matrix2[0]; //v=-ad+bc
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
Erzeugt eine Matrix in der Nähe der Identität, nähert cos(a)=1-a^2/2 & sin(a)=a
*/
complex<double>* set_close_to_identity(complex<double>* matrix){
	double alpha=uniform(0.,delta);
	double u=uniform(-1.,1.);
	double theta=uniform(0., 2.*pi);
	matrix[0]=1.-alpha*alpha/2.+1i*alpha*u;											//cos(a)+isin(a)n_3
	matrix[1]=alpha*(sqrt(1-u*u)*theta+1i*(sqrt(1-u*u)*(1-theta*theta/2.)));		//sin(a)(n_2+in_1)
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
Berechnet Summe der Plaquetts in Vorwärtsrichtung eines Punktes
*/
double calc_plaquette(lattice* lat, int t, int x, int y, int z){
	complex<double>* tmp_matrix = create_identity();
	double result=0;
	int c[4]={t,x,y,z};		//Koordinate des Punkts
	int tc[4]={t,x,y,z};	//Temporäre Koordinate
	int cp1[4];				//Nächster Punkt in Richtung t,x,y,z
	point***** pts = lat->pts;
	//Ermittelt nächste Punkte in Vorwärtsrichtung unter Berücksichtigung der periodischen Randbedingungen
	for(int i=0;i<4;i++){
		if(c[i]==lat->len[i]-1){
			cp1[i]=0;
		}
		else{
			cp1[i]=c[i]+1;
		}
	}
	//Berechnet die Summe der Plaquettes in Vorwärtsrichtung
	for(int mu=0;mu<4;mu++){
		for(int nu=0;nu<4;nu++){
			if(mu<nu){
				//Plaquette in Richtung mu, nu
				tmp_matrix[0]=1;tmp_matrix[1]=0;
				matrix_mul(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],tmp_matrix);		//U_mu(x)
				tc[mu]=cp1[mu];																		//x+a mu
				matrix_mul(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);		//U_nu(x+a mu)
				tc[mu]=c[mu]; tc[nu]=cp1[nu];														//x+a nu
				matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],tmp_matrix);	//U^+_mu(x+a nu)
				tc[nu]=c[nu];																		//x
				matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);	//U^+_nu(x)
				result+=real(2.*(1.-tmp_matrix[0]));
			}
		}
	}
	free(tmp_matrix);
	return result;
}

/*
Berechnet das Staple einer Linkvariable
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
			tmp_matrix[0]=1;tmp_matrix[1]=0;
			tc[mu]=cp1[mu];																		//x+a mu
			matrix_mul(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);		//U_nu(x+a mu)
			tc[mu]=c[mu]; tc[nu]=cp1[nu];														//x+a nu
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],tmp_matrix);	//U^+_mu(x+a nu)
			tc[nu]=c[nu];																		//x
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);	//U^+_nu(x)
			K[0]+=tmp_matrix[0];K[1]+=tmp_matrix[1];
			//Rückwärtsorientierte Schleife
			tmp_matrix[0]=1;tmp_matrix[1]=0;
			tc[mu]=cp1[mu]; tc[nu]=cm1[nu];														//x+a (mu-nu)
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);	//U^+_nu(x+a (mu-nu))
			tc[mu]=c[mu];																		//x-a nu
			matrix_mul_adj(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[mu],tmp_matrix);	//U^+_mu(x-a nu)
			matrix_mul(tmp_matrix,pts[tc[0]][tc[1]][tc[2]][tc[3]]->link[nu],tmp_matrix);		//U(x-a nu)
			K[0]+=tmp_matrix[0]; K[1]+=tmp_matrix[1];
			tc[nu]=c[nu];																		//x
		}
	}
	free(tmp_matrix);
	return K;
}

/*
Berechnet die Änderung in der Wirkung über das Staple K
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
	double S=0;
	//Schleife über alle Punkte in t-Richtung
	for(int t=0;t<len[0];t++){
		//Schleife über alle Punkte in x-Richtung
		for(int x=0;x<len[1];x++){
			//Schleife über alle Punkte in y-Richtung
			for(int y=0;y<len[2];y++){
				//Schleife über alle Punkte in z-Richtung
				for(int z=0;z<len[3];z++){
					S+=calc_plaquette(lat,t,x,y,z);
				}
			}
		}
	}
	return beta/2.*S;
}

/*
Aktualisiert alle Linkvariablen des Gitter N_hit mal
*/
lattice* Metropolis_sweep(lattice* lat){
	complex<double>* tmp_matrix=create_identity();
	complex<double>* K = create_identity();			//Staple zum Punkt in Richtung mu
	int accepted=0; int runs=0;
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
						//Aktualisiert die aktuelle Linkvariable N_hit mal
						for(int i=0;i<N_hit;i++){
							set_close_to_identity(tmp_matrix);
							matrix_mul(pts[t][x][y][z]->link[mu],tmp_matrix,tmp_matrix);
							rescale_su2(tmp_matrix);
							double delta_S=calc_delta_S(tmp_matrix,pts[t][x][y][z]->link[mu],K);
							double r=uniform(0.,1.);
							if((r<exp(-delta_S))){
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
	//Gibt einen Hinweis aus, falls ein Sweep eine Akzeptanz hatte, die nicht zwischen 40% und 60% liegt
	if(((double)accepted/runs<0.4)||((double)accepted/runs>0.6)){
		cout << "Hinweis: Akzeptanz des Sweep betrug " << (double)accepted/runs*100 << "%"<< endl;	
	}
	free(tmp_matrix);
	free(K);
	return lat;
}

/*
Berechnet den Mittelwert einer Observable nach dem Monte-Carlo Verfahren mit Hilfe des Metropolis Sweeps
*/
void Monte_Carlo(){
	clock_t timer=clock();
	double observable[N_cf];
	double estimate=0; double d_estimate=0;
	cout << "Initialisiere Gitter..." << endl;
	lattice* lat=initialize_lattice(N_t,N_s,N_s,N_s);
	cout << "Erfolg" << endl;
	cout << "Thermalisiere Gitter..." << endl;
	for(int i=0;i<N_ini;i++){
		if(i%((int)ceil(N_ini/10.))==0) cout << (double)i*100./N_ini << "%" << endl;
		for(int j=0;j<N_cor+1;j++){
			Metropolis_sweep(lat);
		}
	}
	cout << "Erfolg" << endl;
	cout << "Berechne Werte der Observable..." << endl;
	for(int i=0;i<N_cf;i++){
		for(int j=0;j<N_cor+1;j++){
			Metropolis_sweep(lat);
		}
		observable[i]=calc_action(lat);
		cout << "Wert " << i+1 << ": " << observable[i] << endl;
		estimate+=observable[i];
	}
	estimate=estimate/N_cf;
	for(int i=0;i<N_cf;i++){
		d_estimate+=pow(observable[i]-estimate,2);
	}
	d_estimate=sqrt(d_estimate/(N_cf*(N_cf-1)));
	cout << "Mittelwert: " << estimate << " +/- " << d_estimate << endl;
	free_lattice(lat);
	cout << "Programm nach " << (clock()-timer)/(double) CLOCKS_PER_SEC <<" Sekunden beendet." << endl;
}

//ENDE DER FUNKTIONEN



//BEGINN MAIN

int main(){
	Monte_Carlo();
	return 0;
}
