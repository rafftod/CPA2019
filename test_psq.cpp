#include<iostream>
#include<fstream>
using namespace std;


int main(int argc, char *argv[])
{
	ifstream in ("P00533.fasta.psq", ios::binary); // a changer pour toute la database
	ofstream out ("out.txt");
	if( in.is_open() )
	{
		int8_t x;
		while( in.read((char *) (&x), sizeof(x) )){
			cout << (int) x << endl;
			out << (int) x << endl;// print tout dans un fichier texte pour voir a quoi ca ressemble
		}
		in.close();
		out.close();
int		
	}
	else
		cout << "Impossible d'ouvrir le fichier" << endl;
	return 0;
}
