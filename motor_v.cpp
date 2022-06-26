//Motor v2

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <random>
#include <ctime>
#include <complex>

using std::cout;
using std::endl;
using std::setw;
using cd = std::complex<double>;
double PI = acos(-1);

//criando classe motor v2
class motor {
	
	//especificando os parâmetros de classe
private:
	double l, m, M, r, g, dt, omega, alfa;
	int N;
	//criar construtor da clase motor v2
public:
	motor(double l, double m, double M, double r, double g, double dt, double omega, double alfa, int N) {
		this->g = g;
		this->l= l;
		this->m= m;
		this->M= M;
		this->N= N;
		this->r= r;
		this->dt = dt;
		this->omega = omega;
		this->alfa = alfa;
	}
	
	//criando funções membro
	
	//função para gerar arrays da cinemática
	
	double* tseries(char op) {
		//tempo inicial
		double t = 0.;
		//vetores das posições angulares iniciais
		double* psi = new double[N];
		double* phi = new double[N];
		//vetores das velocidades angulares 
		double* w_phi = new double[N];
		double* w_psi = new double[N];
		double alpha = PI/4.0;
		//condições iniciais
		w_psi[0]= 1.;
		w_phi[0] = 2.0;
		phi[0] = 2.0;
		psi[0] = 1.0;
		
		
		//calculo das séries temporais
		for(int i=0; i<N; i++) {
			//discretização e redução de ordem da EDO	
			w_phi[i+1] =((m*r*cos(omega*t)*omega*omega*l*cos(psi[i])+
								m*r*sin(omega*t)*omega*omega*2*l*sin(psi[i])-m*g*(-sin(alpha)*l*sin(psi[i]))
								+cos(alpha)*l*cos(psi[i])))*sin(alpha)*dt/(m*l*l)+w_phi[i];
								
			w_psi[i+1] = (m*r*cos(omega*t)*omega*omega*l*cos(psi[i])+ m*r*sin(omega*t)*omega*omega*l*sin(psi[i]) 
							- m*g*(-sin(alpha)*l*sin(psi[i]) + cos(alpha)*l*cos(psi[i]))*sin(alpha))*dt/(m*l*l) 
							+ w_psi[i];
			psi[i+1] = 	w_psi[i]*dt + psi[i];		
			phi[i+1] =  w_phi[i]*dt + phi[i];			
			t+= dt;
		}
		
		//instrução switch para retorno de 4 arrays
		
		switch(op) { //op -> pho, psi, wpsi, wphi
			case 'A': {
				return w_phi;
				break;
				
			}
			case 'B': {
				return w_psi;
				break;
			}
			case 'C': {
				return psi;
				break;
				
			}
			
			case 'D': {
				return phi;
				break;
			}
			default:
				cout<<"Escolha incorreta"<<endl;
		}
		
		//liberação de memória
		delete[] phi;
		delete[] psi;
		delete[] w_psi;
		delete[] w_phi;
		return 0;
		
		
	}
	
	//função para imprimir séries temporais
	void print_tsx() {
		std::ofstream output("ts_x.txt");
		char a = 'A';
		double* w_phi = tseries(a);
		double t= 0.0;
		
		for(int i=0;i<N; i++) {
			//cout<<t<<setw(16)<<w_phi[i]<<endl;
			output<<t<<setw(16)<<w_phi[i]<<endl;
			t += dt;
		}
		
	}
	
	void print_tsy() {
		std::ofstream output("ts_y.txt");
		char b = 'B';
		double* w_psi = tseries(b);
		double t= 0.0;
		
		for(int i=0;i<N; i++) {
			//cout<<t<<setw(16)<<w_psi[i]<<endl;
			output<<t<<setw(16)<<w_psi[i]<<endl;
			t += dt;
		}
		
	}
	
	// função que imprime os diagramas de fase
	void print_df(char op){
		switch(op){          // op-> psi, phi, wphi, wpsi
			case 'a':{
				std::ofstream output("psi_vs_wpsi.dat");
				char c = 'C';
				double* psi= tseries(c);
				char b = 'B';
				double* w_psi= tseries(b);
				for(int i=0; i<N; i++){
					output<<psi[i]<<setw(16)<<w_psi[i]<<endl;
					//cout<<psi[i]<<setw(16)<<w_psi[i]<<endl;
					}
					break;
				}
			case 'b':{
				std::ofstream output("phi_vs_wphi.dat");
				char d = 'D';
				double* phi= tseries(d);
				char a = 'A';
				double* w_phi= tseries(a);
				for(int i=0; i<N; i++){
					output<<phi[i]<<setw(16)<<w_phi[i]<<endl;
					//cout<<phi[i]<<setw(16)<<w_phi[i]<<endl;
					}
					break;
				}
			}
		}
	
	
	//transformada de fourier
	void dft_real(char op) {
		double* R = tseries(op);
		std::ofstream output("dft_img.dat");
		double Wmax = PI/dt;
		double dW = 2.*Wmax/(double)(N);
		double Wn = -Wmax;
		
		cd g[N];
		
		double W[N];
		
		cd sum;
		
		for (int i=0; i<N; i++) {
		
			sum = cd(0.0,0.0);
			
			for(int j=0; j<N; j++) {
				
				double tj = dt*j;
				double realPart = cos(Wn*j);
				double imagPart = sin(Wn*tj);
				
				cd Wi(realPart, -imagPart);
				
				sum += R[j]*Wi;
			}
			
			W[i] = Wn;
			
			g[i] = sum/sqrt(double(N)-1.);
			
			Wn += dW;
		
		}
		
		for (int j=0; j<N; j++) {
			output<<W[j]<<setw(16)<<g[j].imag()<<endl;
		}
		
	}
	
};

int main() {
	
	//definindo parametros da classe
	 double l, m, M, r, g, dt, omega, alfa;
	 
	 double T= 50.; //tempo total em segundos
	 dt = 0.01; //passo de tempo
	 omega = 20.; //velocidade angular
	 int N = T/dt;
	 l = 10.0; //comprimento da biela
	 r = 5.0; //comprimento da alavanca
	 m = 5.0; //massa do biela
	 M = 8.0; //massa do pistão
	 g = 9.81; //aceleração gravitacional
	 
	 //criando objeto da classe motor_v
	 
	 motor dyn(l, m, M, r, g, dt, omega, alfa, N);
	 
	 char A, B, C, D;
	 
	 A = 'A';
	 
	 B = 'B';
	 
	 C = 'C';
	 
	 D = 'D';
	 //chamando funções
	 //dyn.print_tsx();
	 
	 //dyn.print_tsy();
	 
	 //dyn.print_df('b'); //-> escolher qual caso quer
	 
	 dyn.dft_real(A);
	 
	 
	
	
	return 0;
}
