#include <iostream>

using namespace std;

#include <complex>
#include <iostream>
#include <valarray>

const double PI = 3.141592653589793238460;

//Declaracion Array de complejos
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;





////Funciones////


//Cool-Tukey FFt (in-place, divide-and-conquer)
//higher memory requirements and redundancy although more intuitive
/*La funcion recive como parametro las samples para la fft, lo que retorna es los
frequency bins de la FFT*/
void fft(CArray &x)
{
    /*Se guarda la referencia a nuestro arreglo de tipo CArray en
    la variable "x", de cierta manera x representa el array dentro de la funcion local*/
    /*obtenemos el size de x y lo guardamos en N, cuano N<=1 significa
    que hemos terminado con la recursion (cada subsuma es de un solo elemento)*/
    const size_t N = x.size();
    if (N<=1) return; //condicion que rompe la recursion


    /*En caso de que N>1 entonces aun no hemos acabado la recursion
    partimos el array en dos subarrays
    uno de index par y otro de index inpar*/
    //devide
    CArray even = x[std::slice(0,N/2,2)];
    CArray odd  = x[std::slice(1,N/2,2)];

    /*recursion: volvemos a llamar a la misma funcion, estavez el parametro
    en vez de ser el array de todas las samples, son los subarrays even y odd*/
    /*En cada iteracion estamos partiendo nuestro array original en 2 arrays
    esto se hace hasta que el size del array es N=1*/
    //conquer
    fft(even);
    fft(odd);
    /*se repiten todas las lineas anteriores hasta que N=1, cuando N=1 termina
    la recursion, la funcion deja de llamarse asi misma y va de regreso ejecutando
    las siguentes lineas*/


    /*se conbinan cada una de las subsumas anteriores la suma*/
    /*esto se hace para cada valor de k(que representa cada frequency bin)*/
    //combine
    for (size_t k=0; k<N/2; k++)
        {
            Complex t = std::polar(1.0, -2*PI*k/N) * odd[k];
            x[k]     = even[k] + t;
            x[k+N/2] = even[k] - t;
        }
}


// inverse fft(in-place)
void ifft(CArray &x)
{
    //conjugate the complex numbers
    x= x.apply(std::conj);

    //forward fft
    fft( x );

    //conjugate the complex numbers again
    x= x.apply(std::conj);

    //scale the numbers
    x /=  x.size();
}


/*other test signals:
const Complex test[]= { 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 };
CArray data(test,8);

const Complex test[]= { 0.0, 1.0, 0.0, -1.0};
CArray data(test,4);
*/

////  Main  ////
int main()
{
    //Test signal
    const Complex test[]= { 0.0, 1.0, 0.0, -1.0};
    CArray data(test,4);

    //other characteristics of the signal
    int N = data.size();                          //Total number of samples
    int nyquistLim = data.size()/2;               //Nyquist Limit
    double Amp[nyquistLim];                       //Amplitude in frequecy spectrum(Averaged final value)
    int samplingFreq= 4;                          //valor dado por el usuario (hz), su unico uso es calcular la resolucion del espectro
    int freqResolution= samplingFreq/data.size(); //especifica el intervalo entre cada freq bin


    //forward fft
    /*recive como parametro las samples almacdnadas en data,
    por que el parametro es pasado por referencia, el resultado es almacenado
    en el mismo arreglo data, esto queire decir que data despues de ejecutar fft(data)
    en lugar de almacenar los samples, ahora almacena las frequency bins resultantes
    esta es una manera de trabajar con arreglos dentro de funciones ya que en c no se puede
    retornar un arreglo entero como tal*/
    fft(data);

    //Average
    /*Obtener la magnitud del numero complejo que representa cada frequecy bin
    una ves obtenido promediarlo para obtener nuestro resultado final
    solo la mitad de las muestras seran tomadas en cuenta por el limite de nyquist*/
    for(int i=0;i<=nyquistLim; i++)
    {
        Amp[i]= std::abs(data[i])*2/data.size();
    }



    //Imprimir FFT (frequency bins)
    std::cout << "FFT (frequency bins complex values): " << std::endl<<"Fk"<<std::endl;

    for (int i=0; i<N; i++)
    {
      std::cout <<"F["<<i<<"]: "<< data[i] << std::endl;
    }
    printf("\n\n");


    //Imprimir freq. Espectrum
    int freq=0;
    std::cout << "Final values[freq spectrum]: " << std::endl;
    std::cout << "Freq.  Amp." << std::endl;
    for (int i=0; i<=nyquistLim; i++)
    {
      std::cout <<freq <<"[Hz]: "<<Amp[i] <<std::endl;
      freq += freqResolution;
    }
    printf("\n\n");


    //inverse fft
    ifft(data);
    std::cout<< std::endl << "ifft"<< std::endl;
    for(int i=0; i<N; i++)
    {
      std::cout << data[i] << std::endl;
    }


    return 0;
}



/* NOTAS
>>Biblioteca estadnar "std"
    En C++, la biblioteca estándar es una colección de clases y funciones, escritas
    en el núcleo del lenguaje. La biblioteca estándar proporciona varios contenedores
    genéricos, funciones para utilizar y manipular esos contenedores, funciones objeto,
    cadenas y flujos genéricos (incluyendo E/S interactiva y de archivos) y soporte para
    la mayoría de las características del lenguaje. La biblioteca estándar de C++ también incorpora la biblioteca estándar de C. Las características de la biblioteca estándar están declaradas en el espacio de nombres std.



>>Clases utilizadas en este programa de la biblioteca std

std::complex
The specializations std::complex<float>, std::complex<double>, and std::complex<long double>
are LiteralTypes for representing and manipulating complex numbers.
Aqui estamos haciendo utilizacion de "templates"


std:valarray<template>


std::slice


std::polar


std::conj


*/



