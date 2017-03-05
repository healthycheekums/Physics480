#include<iostream>

class Complex{
	private:
		//private shiiit
		double real,
				imag;
	public:
		//public shieet -outside the class
		Complex();
		Complex(double);
		Complex(double,double);
		~Complex(){};

		void print();
		Complex operator + (Complex,Complex);
};

Complex::Complex(){
	real= 0.0;
	imag= 0.0;
}
Complex::Complex(double re){
	real=re;
	imag=0.0;
}
Complex::Complex(double re, double im){
	real=re;
	imag=im;
}
void Complex::print(){
	std::cout <<"(" << real << "," << imag << ")" << std::endl;
}
Complex::Complex operator+(Complex a, Complex b){
	return Complex(a.real + b.real, a.imag + b.imag);
}
int main(int argc,char*argv[]){
//declare variables
	double real, imag;

//ask for input
		std::cout<< "Real: ";
		std::cin >> real;
		std::cout <<"Imaginary: ";
		std::cin >> imag;

		Complex number1;
		Complex number2(real);
		Complex number3(real,imag);
		Complex number4 = number2 + number3

		number1.print();
		number2.print();
		number3.print();
		number4.print();

		return 0;
}