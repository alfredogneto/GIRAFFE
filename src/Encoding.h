#pragma once
#include <vector>
#include <string>

#include "Base64.h"

template <typename T>
std::string encodeData(std::vector<T> data);

template <typename T>
std::string encodeData(std::vector<T> data)
{
	int numberBytes = (int)data.size() * 4;	//tamanho do vector * 4 bytes (32 bits)
	char* charNumberBytes = (char*)malloc(8);				//aloca 8 bytes
	char* charData;											//ponteiro não alocado
	char* charFinal = (char*)malloc(8 + data.size() * 4);	//aloca 8 bytes + tamanho do vector*4 bytes
	memset(&charNumberBytes[0], 0, 8);						//memset escreve o valor '0' em cada uma das posições nos 8 bytes
	memset(&charFinal[0], 0, 8 + data.size() * 4);			//memset escreve o valor '0' em cada uma das posições nos 8 bytes + tamanho do vector*4 bytes
	sprintf(charNumberBytes, "%d", numberBytes);
	memcpy(charFinal, charNumberBytes, 8);					//memcpy escreve o string 'charNumberBytes' no inicio do vetor de caracteres 'charFinal'
	for (unsigned int i = 0; i < data.size(); ++i)			//percorre os dados e os coloca no vetor 'charFinal'
	{
		charData = (char*)&data[i];							//inserção do dado proveniente de 'data[i]' (aqui e feito um pointer cast)
		memcpy(&charFinal[8 + i * 4], charData, 4);			//copia para a memória, na posição sequencial, o charData
	}
	std::string encodeFinal = b64encode(charFinal, 8 + data.size() * 4);//realiza encoding com base 64
	free(charFinal);			//desaloca memória
	free(charNumberBytes);		//desaloca memória
	return encodeFinal;
}
