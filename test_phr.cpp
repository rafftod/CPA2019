#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

int main(int argc, char *argv[])
{
	ifstream in ("uniprot_sprot.fasta.phr", ios::binary);
	ofstream out ("out_phr.txt");
	
	if( in.is_open() )
	{
		uint8_t x;
		
		while( in.read((char *) (&x), sizeof(x) )){ //lecture byte par byte (sizeof(x) = 1 )
			cout << hex << (int32_t) x << endl;
		}
	}
	else
		cout << "Impossible d'ouvrir le fichier" << endl;
		
	in.close();
	out.close();
	cout << endl;
	
	return 0;
	
}
