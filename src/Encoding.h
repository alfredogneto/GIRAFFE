#pragma once
#include <vector>
#include <string>

#include "Base64.h"

template <typename T>
std::string encodeData(std::vector<T> data);

template <typename T>
std::string encodeData(std::vector<T> data)
{
	int numberBytes = (int)data.size() * 4;	//size of the vector * 4 bytes (32 bits)
	char* charNumberBytes = (char*)malloc(8);				//allocates 8 bytes
	char* charData;											//pointer not allocated
	char* charFinal = (char*)malloc(8 + data.size() * 4);	//allocates 8 bytes + size of the vector*4 bytes
	memset(&charNumberBytes[0], 0, 8);						//memset writes the value '0' in each a of the positions in the 8 bytes
	memset(&charFinal[0], 0, 8 + data.size() * 4);			//memset writes the value '0' in each a of the positions in the 8 bytes + size of the vector*4 bytes
	sprintf(charNumberBytes, "%d", numberBytes);
	memcpy(charFinal, charNumberBytes, 8);					//memcpy writes the string 'charNumberBytes' in the start of the vector of caracteres 'charFinal'
	for (unsigned int i = 0; i < data.size(); ++i)			//iterates over the data and the coloca in the vector 'charFinal'
	{
		charData = (char*)&data[i];							//insertion of the data proveniente of 'data[i]' (here and done a pointer cast)
		memcpy(&charFinal[8 + i * 4], charData, 4);			//copies for the memory, in the position sequencial, the charData
	}
	std::string encodeFinal = b64encode(charFinal, 8 + data.size() * 4);//performs encoding with base 64
	free(charFinal);			//deallocates memory
	free(charNumberBytes);		//deallocates memory
	return encodeFinal;
}
