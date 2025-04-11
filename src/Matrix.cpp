#include "Matrix.h"
#include <mkl.h>

#define PI 3.1415926535897932384626433832795

//Construtor Padrão
Matrix::Matrix(void)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = 1;
	m_columns = 1;
	//Inicializa a matriz como null e tenta aloca-la
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Construtor de matriz quadrada
Matrix::Matrix(long lines)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = lines;
	m_columns = 1;
	//Inicializa a matriz como null e tenta aloca-la
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Construtor Parametrico
Matrix::Matrix(long lines, long columns)
{
	m_lines_deleted = true;
	m_alloced_lines = 0;
	m_lines = lines;
	m_columns = columns;
	//Inicializa a matriz como null e tenta aloca-la
	m_matrix = NULL;
	if(!alloc())
		printf("Nao foi possivel alocar matriz! \n");
}
//Construtor de cópia
Matrix::Matrix(const Matrix &copied)
{
	//Verifica dimensões da matriz - se necessario, faz re-alocação
	//if (copied.m_alloced_lines != m_alloced_lines)
	//{
		m_alloced_lines = 0;
		m_lines = copied.m_lines;
		m_columns = copied.m_columns;
		m_lines_deleted = true;
		//Inicializa a matriz como null e tenta aloca-la
		m_matrix = NULL;
		if(!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	//}
	//Copia os valores na nova matriz
	for(long i=0; i < copied.m_alloced_lines; i++)
		m_matrix[i] = copied.m_matrix[i];	
}
//Destrutor Padrão
Matrix::~Matrix(void)
{
	flush();
}

//Retorna o numero de linhas da matriz
long Matrix::getLines()
{
	return this->m_lines;
}
//Retorna o numero de colunas da matriz 
long Matrix::getColumns()
{
	return this->m_columns;
}
//Define o numero de linhas da matriz
void Matrix::setLines(long value) 
{
	m_lines = value;
}
//Define o numero de colunas da matriz
void Matrix::setColumns(long value) 
{
	m_columns = value;
}

//Imprime a matriz na tela
void Matrix::print()
{
	printf("\n");
	for(long i=0; i < this->getLines(); i++)
	{
		printf("|");
		for(long j=0; j < this->getColumns(); j++)
		{
			if(m_matrix[i+j*m_lines] >= 0.0)
				printf(" %.4e ",m_matrix[i+j*m_lines]);
			else
				printf(" %.4e ",m_matrix[i+j*m_lines]);
		}
		printf("|\n");
	}
	printf("\n");
}
//Imprime a matriz em um arquivo de texto, cujo nome esta no char s
void Matrix::fprint(char* s)
{
	FILE *file1 = fopen(s,"w");
	fprintf(file1,"\n");
	for(long i=0; i < this->getLines(); i++)
	{
		for(long j=0; j < this->getColumns(); j++)
		{
			if(m_matrix[i+j*m_lines] >= 0.0)
				fprintf(file1," %.14e\t",m_matrix[i+j*m_lines]);
			else
				fprintf(file1," %.14e\t",m_matrix[i+j*m_lines]);
		}
		fprintf(file1,"\n");
	}
	fprintf(file1,"\n");
	fclose(file1);
}
//Aloca a matriz
bool Matrix::alloc()
{
	flush();
	//Tenta alocar a matriz
	if(!(m_matrix = new double[m_lines*m_columns]))
			//Falha
			return false;
	//Caso consiga
	else
	{
		m_alloced_lines = m_lines*m_columns;
		this->m_lines_deleted = false;
		
		for(long i=0; i < m_alloced_lines; i++)
				m_matrix[i] = 0.0;
		return true;
	}	
}
//Libera a memória ocupada pela matriz
bool Matrix::flush()
{
	//Percorre as linhas
	//Se ha linhas para desalocar
	if (this->m_lines_deleted == false)
	{
		if (m_matrix)
			delete[]m_matrix;
		//Desalocada com sucesso
		m_lines_deleted=true;	//para evitar desalocação novamente
		m_matrix = NULL;
		m_alloced_lines = 0;
	}
	return true;
}
//Zera a matriz, mantando as dimensões atuais
void Matrix::clear()
{
	for(int i=0; i < m_alloced_lines; i++)
			m_matrix[i] = 0.0;
}

//Operador Soma
Matrix operator + (Matrix &matrix1, Matrix &matrix2)
{
	//Verifica se as dimensões das matrizes são compativeis
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem de erro - Dimensões Incompativeis
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Retorna vazio
		return NULL;
	}
	//Caso as dimensões sejam compativeis
	else
	{
		//Cria uma matriz de retorno
		Matrix return_matrix(matrix1.getLines(),matrix1.getColumns());
		
		//Subtrai elemento a elemento da matriz
		for(long i=0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = matrix1.m_matrix[i] + matrix2.m_matrix[i];
		
		//Retorna a matriz subtração
		return return_matrix;
	}
}
//Operador Subtração
Matrix operator - (Matrix &matrix1, Matrix &matrix2)
{
	//Verifica se as dimensões das matrizes são compativeis
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
	{
		//Mensagem de erro - Dimensões Incompativeis
		printf("Matrizes devem possuir a mesma dimensao! \n");
		//Retorna vazio
		return NULL;
	}
	//Caso as dimensões sejam compativeis
	else
	{
		//Cria uma matriz de retorno
		Matrix return_matrix(matrix1.getLines(),matrix1.getColumns());
		
		//Subtrai elemento a elemento da matriz
		for(long i=0; i < matrix1.m_alloced_lines; i++)
			return_matrix.m_matrix[i] = matrix1.m_matrix[i] - matrix2.m_matrix[i];
		
		//Retorna a matriz subtração
		return return_matrix;
	}
}
//Operador Multiplicacao de matrizes
Matrix operator * (Matrix &matrix1, Matrix &matrix2)
{
	//Verificação da possibilidade de multiplicação
	if (matrix1.m_columns != matrix2.m_lines)
	{
		//Produto escalar
		if (matrix1.m_lines == matrix2.m_lines)
		{
			Matrix ret(1,1);
			for (int i = 0; i < matrix1.m_lines; i++)
				ret(0,0) += matrix1(i, 0)*matrix2(i, 0);
			return ret;
		}
		else
		{
			printf("Nao e possivel multiplicar as matrizes. Dimensoes incompativeis!");
			return NULL;
		}
	}
	else
	{
		
		Matrix return_m(matrix1.m_lines,matrix2.m_columns);
		//Se a segunda matriz for um vetor
		if (matrix2.m_columns == 1)
			cblas_dgemv(CblasColMajor,CblasNoTrans,matrix1.m_lines,matrix1.m_columns,1.0,matrix1.m_matrix,matrix1.m_lines,matrix2.m_matrix,1,0.0,return_m.m_matrix,1);
		//Se não
		else
			cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,matrix1.m_lines,matrix2.m_columns,matrix2.m_lines,1.0,matrix1.m_matrix,matrix1.m_lines,matrix2.m_matrix,matrix2.m_lines,1,return_m.m_matrix,return_m.m_lines);
		//return_m.print();
		return return_m;
	}
}
//Operador Multiplicacao por escalar
Matrix operator * (double escalar, Matrix &matrix1)
{
	//Cria uma matriz de retorno
	Matrix return_matrix(matrix1.getLines(),matrix1.getColumns());

	//Faz o produto de cada elemento da matriz pelo escalar
	for(long i=0; i < return_matrix.m_alloced_lines; i++)
		return_matrix.m_matrix[i] = matrix1.m_matrix[i]*escalar;
	
	//Retorna a matriz produto por escalar
	return return_matrix;

	
}
//Operador Multiplicacao por escalar
Matrix operator * (Matrix &matrix1, double escalar)
{
	//Retorna a função anterior, ja que esta operação e comutativa
	return escalar*matrix1;
}
//Verificação de igualdade
bool operator == (Matrix &matrix1, Matrix &matrix2)
{
	//Verifica a igualdade de dimensões
	if((matrix1.getLines() != matrix2.getLines()) || (matrix1.getColumns() != matrix2.getColumns()))
		return false;
	else
	{
		//Varre a matriz e verifica a igualdade de elementos
		for(long i=0; i < matrix1.m_alloced_lines; i++)
				if(matrix1.m_matrix[i] != matrix2.m_matrix[i]) return false;
		//Matrizes iguais
		return true;
	}
}
//Verificação de inegualdade
bool operator != (Matrix &matrix1, Matrix &matrix2)
{
	return !(matrix1 == matrix2);
}
//Operador de Atribuição	
Matrix &Matrix::operator = (Matrix const &matrix1)
{
	//Verifica dimensões da matriz - se necessario, faz re-alocação
	if (matrix1.m_alloced_lines != m_alloced_lines)
	{
		this->flush();
		m_alloced_lines = 0;
		m_lines = matrix1.m_lines;
		m_columns = matrix1.m_columns;
		m_lines_deleted = matrix1.m_lines_deleted;
		//Inicializa a matriz como null e tenta aloca-la
		m_matrix = NULL;
		if(!alloc())
			printf("Nao foi possivel alocar matriz! \n");
	}
	//Copia os valores na nova matriz
	for(long i=0; i < matrix1.m_alloced_lines; i++)
		m_matrix[i] = matrix1.m_matrix[i];

	//Retorna esta matriz
	return *this;
}
//Retorno do valor na posição especificada
double &Matrix::operator() (long line, long column)
{
	
	//Verifica se a posição e valida
	if(line > this->m_lines-1 || column > this->m_columns-1 || line < 0 || column < 0)
	{
		printf("Not valid position accessed in matrix! (%d,%d)\n",line,column);
		double* ret = new double[1];
		ret[0] = 0;
		return ret[0];
	}
	else
		//Retorna o valor na posição desejada
		return this->m_matrix[line + column*m_lines];
}
//Operador produto escalar entre dois vetores
double dot(Matrix &matrix1, Matrix &matrix2)
{
	//Verificação da possibilidade do produto
	if (matrix1.m_lines != matrix2.m_lines)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		double return_value = 0.0;
		for (int i = 0; i < matrix1.m_alloced_lines; i++)
			return_value+=matrix1.m_matrix[i]*matrix2.m_matrix[i];
		return return_value;
	}
}
//Operador produto vetorial entre dois vetores
Matrix cross(Matrix &matrix1, Matrix &matrix2)
{
	//Verificação da possibilidade do produto
	if (matrix1.m_columns != 1 || matrix2.m_columns != 1 || matrix1.m_lines != 3 || matrix2.m_lines != 3)
	{
		printf("Nao e possivel calcular o produto vetorial. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		Matrix return_m(3);
		return_m(0,0) = matrix1(1,0)*matrix2(2,0) - matrix1(2,0)*matrix2(1,0);
		return_m(1,0) = matrix1(2,0)*matrix2(0,0) - matrix1(0,0)*matrix2(2,0);
		return_m(2,0) = matrix1(0,0)*matrix2(1,0) - matrix1(1,0)*matrix2(0,0);
		return return_m;
	}
}
//Inverte uma matriz 2x2
Matrix invert2x2(Matrix &matrix)
{
	//Verificação da possibilidade de inversão
	if (matrix.getColumns() != 2 || matrix.getLines() != 2)
	{
		printf("Nao e possivel inverter a matriz. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		double a = matrix(0,0);
		double b = matrix(0,1);
		double c = matrix(1,0);
		double d = matrix(1,1);
		double determ = a*d - b*c;
		Matrix return_m(2,2);
		return_m(0,0) = d/determ;
		return_m(0,1) = -b/determ;
		return_m(1,0) = -c/determ;
		return_m(1,1) = a/determ;
	
		return return_m;
	}
}

//Inverte uma matriz 3x3
Matrix invert3x3(Matrix &matrix)
{
	//Verificação da possibilidade de inversão
	if (matrix.getColumns() != 3 || matrix.getLines() != 3)
	{
		printf("Nao e possivel inverter a matriz. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		double determinant = matrix(0,0)*matrix(1,1)*matrix(2,2) + matrix(1,0)*matrix(2,1)*matrix(0,2) + matrix(2,0)*matrix(0,1)*matrix(1,2) - 
							 matrix(0,0)*matrix(2,1)*matrix(1,2) - matrix(2,0)*matrix(1,1)*matrix(0,2) - matrix(1,0)*matrix(0,1)*matrix(2,2);
		
		Matrix return_m(3,3);
		return_m(0,0) = 1.0/determinant*( matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1));
		return_m(0,1) = 1.0/determinant*( matrix(0,2)*matrix(2,1) - matrix(0,1)*matrix(2,2));
		return_m(0,2) = 1.0/determinant*( matrix(0,1)*matrix(1,2) - matrix(0,2)*matrix(1,1));

		return_m(1,0) = 1.0/determinant*( matrix(1,2)*matrix(2,0) - matrix(1,0)*matrix(2,2));
		return_m(1,1) = 1.0/determinant*( matrix(0,0)*matrix(2,2) - matrix(0,2)*matrix(2,0));
		return_m(1,2) = 1.0/determinant*( matrix(0,2)*matrix(1,0) - matrix(0,0)*matrix(1,2));

		return_m(2,0) = 1.0/determinant*( matrix(1,0)*matrix(2,1) - matrix(1,1)*matrix(2,0));
		return_m(2,1) = 1.0/determinant*( matrix(0,1)*matrix(2,0) - matrix(0,0)*matrix(2,1));
		return_m(2,2) = 1.0/determinant*( matrix(0,0)*matrix(1,1) - matrix(0,1)*matrix(1,0));

		return return_m;
	}
}

//Inverte uma matriz 4x4
Matrix invert4x4(Matrix &matrix)
{
	double v[198];
	Matrix invM(4, 4);
	v[1] = matrix(0,0);
	v[2] = matrix(0,1);
	v[3] = matrix(0,2);
	v[4] = matrix(0,3);
	v[5] = matrix(1,0);
	v[6] = matrix(1,1);
	v[85] = -(v[2] * v[5]) + v[1] * v[6];
	v[7] = matrix(1,2);
	v[88] = v[3] * v[6] - v[2] * v[7];
	v[79] = v[3] * v[5] - v[1] * v[7];
	v[60] = v[4] * v[7];
	v[8] = matrix(1,3);
	v[84] = -(v[4] * v[6]) + v[2] * v[8];
	v[78] = -(v[4] * v[5]) + v[1] * v[8];
	v[59] = v[3] * v[8];
	v[77] = -v[59] + v[60];
	v[9] = matrix(2,0);
	v[10] = matrix(2,1);
	v[86] = v[10] * v[7];
	v[83] = -(v[1] * v[10]) + v[2] * v[9];
	v[82] = v[10] * v[5] - v[6] * v[9];
	v[80] = v[10] * v[8];
	v[11] = matrix(2,2);
	v[87] = v[11] * v[6];
	v[76] = v[1] * v[11] - v[3] * v[9];
	v[74] = -(v[11] * v[5]) + v[7] * v[9];
	v[64] = v[11] * v[8];
	v[62] = v[11] * v[4];
	v[12] = matrix(2,3);
	v[81] = v[12] * v[6];
	v[75] = -(v[1] * v[12]) + v[4] * v[9];
	v[73] = v[12] * v[5] - v[8] * v[9];
	v[63] = v[12] * v[7];
	v[61] = v[12] * v[3];
	v[72] = v[10] * (-v[59] + v[60]) + v[6] * (v[61] - v[62]) + v[2] * (-v[63] + v[64]);
	v[13] = matrix(3,0);
	v[14] = matrix(3,1);
	v[15] = matrix(3,2);
	v[66] = v[15] * v[4];
	v[65] = v[15] * v[2];
	v[16] = matrix(3,3);
	v[68] = v[16] * v[3];
	v[67] = v[16] * v[2];
	v[71] = v[14] * (v[59] - v[60]) + v[6] * (v[66] - v[68]) + v[67] * v[7] - v[65] * v[8];
	v[70] = v[14] * (-v[61] + v[62]) + v[12] * v[65] - v[11] * v[67] + v[10] * (-v[66] + v[68]);
	v[69] = v[14] * (v[63] - v[64]) + v[15] * (v[80] - v[81]) + v[16] * (-v[86] + v[87]);
	v[17] = 1e0 / (v[1] * v[69] + v[5] * v[70] + v[13] * v[72] + v[71] * v[9]);
	invM(0,0) = v[17] * v[69];
	invM(0,1) = v[17] * v[70];
	invM(0,2) = v[17] * v[71];
	invM(0,3) = v[17] * v[72];
	invM(1,0) = v[17] * (v[13] * (-v[63] + v[64]) + v[15] * v[73] + v[16] * v[74]);
	invM(1,1) = v[17] * (v[13] * (v[61] - v[62]) + v[15] * v[75] + v[16] * v[76]);
	invM(1,2) = v[17] * (v[13] * v[77] + v[15] * v[78] + v[16] * v[79]);
	invM(1,3) = v[17] * (-(v[11] * v[78]) - v[12] * v[79] - v[77] * v[9]);
	invM(2,0) = v[17] * (-(v[14] * v[73]) + v[13] * (-v[80] + v[81]) + v[16] * v[82]);
	invM(2,1) = v[17] * (v[13] * (-(v[12] * v[2]) + v[10] * v[4]) - v[14] * v[75] + v[16] * v[83]);
	invM(2,2) = v[17] * (-(v[14] * v[78]) + v[13] * v[84] + v[16] * v[85]);
	invM(2,3) = v[17] * (v[10] * v[78] - v[12] * v[85] - v[84] * v[9]);
	invM(3,0) = v[17] * (-(v[14] * v[74]) - v[15] * v[82] + v[13] * (v[86] - v[87]));
	invM(3,1) = v[17] * (v[13] * (v[11] * v[2] - v[10] * v[3]) - v[14] * v[76] - v[15] * v[83]);
	invM(3,2) = v[17] * (-(v[14] * v[79]) - v[15] * v[85] + v[13] * v[88]);
	invM(3,3) = v[17] * (v[10] * v[79] + v[11] * v[85] - v[88] * v[9]);

	return invM;
}

//Inverte uma matriz 5x5
Matrix invert5x5(Matrix &matrix)
{
	double v[481];
	Matrix invM(5, 5);

	v[1] = matrix(0,0);
	v[2] = matrix(0,1);
	v[3] = matrix(0,2);
	v[4] = matrix(0,3);
	v[5] = matrix(0,4);
	v[6] = matrix(1,0);
	v[322] = v[3] * v[6];
	v[7] = matrix(1,1);
	v[306] = v[4] * v[7];
	v[8] = matrix(1,2);
	v[323] = v[1] * v[8];
	v[313] = v[5] * v[8];
	v[301] = v[4] * v[8];
	v[9] = matrix(1,3);
	v[307] = v[2] * v[9];
	v[302] = v[3] * v[9];
	v[247] = v[5] * v[9];
	v[10] = matrix(1,4);
	v[314] = v[10] * v[3];
	v[246] = v[10] * v[4];
	v[11] = matrix(2,0);
	v[12] = matrix(2,1);
	v[13] = matrix(2,2);
	v[321] = v[1] * v[13];
	v[320] = v[13] * v[6];
	v[14] = matrix(2,3);
	v[251] = v[10] * v[14];
	v[249] = v[14] * v[5];
	v[117] = v[2] * v[251];
	v[116] = -(v[249] * v[7]);
	v[15] = matrix(2,4);
	v[250] = v[15] * v[9];
	v[248] = v[15] * v[4];
	v[119] = -(v[2] * v[250]);
	v[118] = v[248] * v[7];
	v[339] = v[116] + v[117] + v[118] + v[119] + v[12] * (-v[246] + v[247]);
	v[16] = matrix(3,0);
	v[315] = v[15] * v[16];
	v[312] = v[13] * v[16];
	v[17] = matrix(3,1);
	v[300] = v[14] * v[17];
	v[299] = v[13] * v[17];
	v[89] = -(v[17] * v[248]);
	v[88] = v[17] * v[249];
	v[76] = v[17] * v[250];
	v[75] = -(v[17] * v[251]);
	v[18] = matrix(3,2);
	v[318] = v[1] * v[18];
	v[317] = v[18] * v[6];
	v[316] = v[11] * v[18];
	v[304] = v[14] * v[18];
	v[303] = v[12] * v[18];
	v[19] = matrix(3,3);
	v[305] = v[19] * v[7];
	v[274] = v[12] * v[19];
	v[263] = v[15] * v[19];
	v[253] = v[19] * v[5];
	v[252] = v[19] * v[2];
	v[104] = -(v[10] * v[252]);
	v[103] = v[253] * v[7];
	v[90] = -(v[12] * v[253]);
	v[77] = v[10] * v[274];
	v[20] = matrix(3,4);
	v[319] = v[11] * v[20];
	v[262] = v[14] * v[20];
	v[258] = v[12] * v[20];
	v[255] = v[20] * v[4];
	v[254] = v[2] * v[20];
	v[106] = v[254] * v[9];
	v[105] = -(v[255] * v[7]);
	v[334] = v[103] + v[104] + v[105] + v[106] + v[17] * (v[246] - v[247]);
	v[92] = v[12] * v[255];
	v[333] = v[15] * v[252] - v[14] * v[254] + v[88] + v[89] + v[90] + v[92];
	v[79] = -(v[258] * v[9]);
	v[328] = (v[262] - v[263])*v[7] + v[75] + v[76] + v[77] + v[79];
	v[199] = v[1] * v[328] + v[11] * v[334] + v[16] * v[339] + v[333] * v[6];
	v[21] = matrix(4,0);
	v[265] = v[18] * v[21];
	v[264] = v[13] * v[21];
	v[261] = v[21] * v[8];
	v[260] = v[21] * v[3];
	v[259] = v[17] * v[21];
	v[257] = v[21] * v[7];
	v[256] = v[2] * v[21];
	v[217] = v[256] * v[8];
	v[216] = v[257] * v[3];
	v[204] = v[12] * v[265];
	v[202] = v[13] * v[259];
	v[155] = v[260] * v[9];
	v[154] = v[261] * v[4];
	v[151] = v[21] * v[246];
	v[150] = v[21] * v[247];
	v[135] = v[21] * v[262];
	v[133] = v[21] * v[263];
	v[132] = v[19] * v[264];
	v[130] = v[14] * v[265];
	v[22] = matrix(4,1);
	v[273] = v[18] * v[22];
	v[272] = v[13] * v[22];
	v[271] = v[22] * v[8];
	v[270] = v[22] * v[3];
	v[269] = v[16] * v[22];
	v[268] = v[11] * v[22];
	v[267] = v[22] * v[6];
	v[266] = v[1] * v[22];
	v[187] = v[10] * v[266];
	v[186] = v[267] * v[5];
	v[175] = v[20] * v[268];
	v[173] = v[15] * v[269];
	v[113] = v[270] * v[9];
	v[112] = v[271] * v[4];
	v[111] = v[10] * v[270];
	v[110] = v[271] * v[5];
	v[335] = -v[110] + v[111];
	v[109] = v[22] * v[246];
	v[108] = v[22] * v[247];
	v[87] = v[22] * v[262];
	v[86] = v[20] * v[272];
	v[85] = v[22] * v[263];
	v[84] = v[19] * v[272];
	v[83] = v[15] * v[273];
	v[329] = -v[83] + v[86];
	v[82] = v[14] * v[273];
	v[23] = matrix(4,2);
	v[281] = v[16] * v[23];
	v[280] = v[11] * v[23];
	v[279] = v[23] * v[6];
	v[278] = v[1] * v[23];
	v[277] = v[17] * v[23];
	v[276] = v[2] * v[23];
	v[275] = v[23] * v[7];
	v[221] = v[1] * v[275];
	v[220] = v[276] * v[6];
	v[369] = v[216] - v[217] + v[220] - v[221] - v[267] * v[3] + v[266] * v[8];
	v[363] = -v[216] + v[217] - v[220] + v[221];
	v[207] = v[11] * v[277];
	v[206] = v[12] * v[281];
	v[366] = v[202] - v[204] + v[206] - v[207] + v[18] * v[268] - v[13] * v[269];
	v[358] = -v[202] + v[204] - v[206] + v[207];
	v[161] = v[278] * v[9];
	v[160] = v[279] * v[4];
	v[344] = -v[154] + v[155] + v[160] - v[161];
	v[138] = v[19] * v[280];
	v[136] = v[14] * v[281];
	v[340] = -v[130] + v[132] + v[136] - v[138];
	v[24] = matrix(4,3);
	v[289] = v[17] * v[24];
	v[288] = v[18] * v[24];
	v[287] = v[24] * v[7];
	v[286] = v[2] * v[24];
	v[285] = v[16] * v[24];
	v[284] = v[11] * v[24];
	v[283] = v[24] * v[6];
	v[282] = v[1] * v[24];
	v[371] = v[154] - v[155] - v[160] + v[161] + v[283] * v[3] - v[282] * v[8];
	v[368] = v[130] - v[132] - v[136] + v[138] - v[18] * v[284] + v[13] * v[285];
	v[353] = v[2] * v[283] - v[267] * v[4] - v[282] * v[7] + v[266] * v[9];
	v[348] = v[19] * v[268] - v[14] * v[269] - v[17] * v[284] + v[12] * v[285];
	v[125] = v[286] * v[8];
	v[124] = v[287] * v[3];
	v[370] = -v[112] + v[113] - v[124] + v[125] + v[275] * v[4] - v[276] * v[9];
	v[338] = v[112] - v[113] + v[124] - v[125];
	v[123] = v[10] * v[286];
	v[122] = v[287] * v[5];
	v[337] = v[108] - v[109] - v[122] + v[123];
	v[121] = v[24] * v[314];
	v[120] = v[24] * v[313];
	v[336] = v[120] - v[121];
	v[99] = v[13] * v[20] * v[24];
	v[98] = v[24] * v[258];
	v[97] = v[15] * v[288];
	v[332] = v[97] - v[99];
	v[96] = v[12] * v[288];
	v[95] = v[15] * v[289];
	v[331] = v[85] - v[87] - v[95] + v[98];
	v[94] = v[13] * v[289];
	v[367] = -(v[23] * v[274]) + v[14] * v[277] - v[82] + v[84] - v[94] + v[96];
	v[330] = v[82] - v[84] + v[94] - v[96];
	v[25] = matrix(4,4);
	v[298] = v[16] * v[25];
	v[297] = v[18] * v[25];
	v[296] = v[19] * v[25];
	v[295] = v[25] * v[6];
	v[294] = v[1] * v[25];
	v[293] = v[12] * v[25];
	v[292] = v[17] * v[25];
	v[291] = v[2] * v[25];
	v[290] = v[25] * v[7];
	v[364] = v[110] - v[111] + v[10] * v[276] + v[290] * v[3] - v[275] * v[5] - v[291] * v[8];
	v[359] = v[23] * v[258] - v[15] * v[277] + v[13] * v[292] - v[18] * v[293] + v[83] - v[86];
	v[197] = v[1] * v[290];
	v[196] = v[291] * v[6];
	v[362] = -v[186] + v[187] + v[196] - v[197] - v[10] * v[256] + v[257] * v[5];
	v[355] = v[186] - v[187] - v[196] + v[197];
	v[354] = -v[108] + v[109] + v[122] - v[123] - v[290] * v[4] + v[291] * v[9];
	v[351] = v[25] * v[274] - v[14] * v[292] - v[85] + v[87] + v[95] - v[98];
	v[181] = v[11] * v[292];
	v[180] = v[16] * v[293];
	v[357] = -v[173] + v[175] + v[180] - v[181] - v[21] * v[258] + v[15] * v[259];
	v[350] = v[173] - v[175] - v[180] + v[181];
	v[347] = v[10] * (-v[260] + v[278]) + v[295] * v[3] + (v[261] - v[279])*v[5] - v[294] * v[8];
	v[165] = v[294] * v[9];
	v[164] = v[295] * v[4];
	v[352] = v[150] - v[151] + v[164] - v[165] + v[10] * v[282] - v[283] * v[5];
	v[346] = -v[150] + v[151] - v[164] + v[165];
	v[345] = -v[120] + v[121] + v[23] * (-v[246] + v[247]) + v[25] * (v[301] - v[302]);
	v[146] = v[11] * v[296];
	v[343] = v[23] * (-v[262] + v[263]) - v[13] * v[296] + v[14] * v[297] - v[97] + v[99];
	v[143] = v[14] * v[298];
	v[349] = v[133] - v[135] + v[143] - v[146] + v[20] * v[284] - v[15] * v[285];
	v[342] = -v[133] + v[135] - v[143] + v[146];
	v[341] = v[20] * (-v[264] + v[280]) + v[15] * (v[265] - v[281]) - v[11] * v[297] + v[13] * v[298];
	v[28] = v[300] * v[8];
	v[29] = v[3] * v[300];
	v[36] = v[18] * v[306];
	v[37] = v[18] * v[307];
	v[38] = v[274] * v[8];
	v[39] = v[274] * v[3];
	v[308] = v[28] + v[13] * v[305] - v[38] - v[304] * v[7] + (-v[299] + v[303])*v[9];
	v[309] = -(v[13] * v[252]) - v[29] + v[2] * v[304] + v[39] + (v[299] - v[303])*v[4];
	v[310] = v[17] * (-v[301] + v[302]) - v[3] * v[305] + v[36] - v[37] + v[252] * v[8];
	v[46] = v[13] * v[306];
	v[47] = v[13] * v[307];
	v[311] = v[12] * (v[301] - v[302]) - v[46] + v[47] + v[14] * (v[3] * v[7] - v[2] * v[8]);
	v[243] = v[1] * v[308] + v[11] * v[310] + v[16] * v[311] + v[309] * v[6];
	v[127] = -(v[10] * v[309]) - v[15] * v[310] - v[20] * v[311] - v[308] * v[5];
	v[50] = v[312] * v[5];
	v[51] = v[10] * v[312];
	v[54] = v[315] * v[8];
	v[55] = v[3] * v[315];
	v[56] = v[316] * v[5];
	v[360] = -v[50] + v[56];
	v[57] = v[10] * v[316];
	v[356] = v[51] - v[57];
	v[58] = v[317] * v[5];
	v[59] = v[10] * v[318];
	v[64] = v[320] * v[5];
	v[65] = v[10] * v[321];
	v[66] = v[15] * v[322];
	v[67] = v[15] * v[323];
	v[365] = -v[66] + v[67];
	v[324] = v[11] * (-v[313] + v[314]) + v[365] + v[64] - v[65];
	v[68] = v[319] * v[8];
	v[69] = v[3] * v[319];
	v[325] = v[15] * v[317] - v[20] * v[320] + v[356] - v[54] + v[68];
	v[326] = -(v[15] * v[318]) + v[20] * v[321] + v[360] + v[55] - v[69];
	v[72] = v[20] * v[322];
	v[73] = v[20] * v[323];
	v[361] = v[72] - v[73];
	v[327] = v[16] * (v[313] - v[314]) + v[361] - v[58] + v[59];
	v[225] = v[17] * v[324] + v[2] * v[325] + v[12] * v[327] + v[326] * v[7];
	v[169] = -(v[19] * v[324]) - v[14] * v[327] - v[325] * v[4] - v[326] * v[9];
	v[74] = 1e0 / (v[127] * v[21] - v[169] * v[22] + v[199] * v[23] - v[225] * v[24] + v[243] * v[25]);
	invM(0,0) = v[74] * (v[25] * v[308] + v[23] * v[328] + v[10] * v[330] + v[332] * v[7] + v[331] * v[8] + v[329] * v[9]);
	invM(0,1) = (v[25] * v[309] - v[3] * v[331] - v[2] * v[332] + v[23] * v[333] - v[329] * v[4] - v[330] * v[5])*v[74];
	invM(0,2) = (v[25] * v[310] + v[23] * v[334] + v[19] * v[335] + v[17] * v[336] + v[18] * v[337] + v[20] * v[338])*v[74];
	invM(0,3) = (v[25] * v[311] - v[14] * v[335] - v[12] * v[336] - v[13] * v[337] - v[15] * v[338] + v[23] * v[339])*v[74];
	invM(0,4) = v[127] * v[74];
	invM(1,0) = v[74] * (v[10] * v[340] + v[343] * v[6] + v[24] * (-v[51] + v[54] + v[57] - v[68]) + v[342] * v[8] + v[341] * v[9]
		);
	invM(1,1) = (-(v[3] * v[342]) - v[1] * v[343] - v[341] * v[4] - v[340] * v[5] + v[24] * (v[50] - v[55] - v[56] + v[69])
		)*v[74];
	invM(1,2) = (v[20] * v[344] + v[16] * v[345] + v[18] * v[346] + v[19] * v[347] + v[24] * (v[58] - v[59] - v[72] + v[73])
		)*v[74];
	invM(1,3) = (-(v[15] * v[344]) - v[11] * v[345] - v[13] * v[346] - v[14] * v[347] + v[24] * (-v[64] + v[65] + v[66] - v[67])
		)*v[74];
	invM(1,4) = -(v[169] * v[74]);
	invM(2,0) = v[74] * (v[10] * v[348] + v[351] * v[6] + v[349] * v[7] + v[21] * (-v[75] - v[76] - v[77] - v[79]) + v[350] * v[9]
		);
	invM(2,1) = v[74] * (-(v[2] * v[349]) - v[1] * v[351] - v[350] * v[4] - v[348] * v[5] + v[21] * (-v[88] - v[89] - v[90]
		- v[92]));
	invM(2,2) = ((-v[103] - v[104] - v[105] - v[106])*v[21] + v[17] * v[352] + v[20] * v[353] + v[16] * v[354]
		+ v[19] * v[355])*v[74];
	invM(2,3) = ((-v[116] - v[117] - v[118] - v[119])*v[21] - v[12] * v[352] - v[15] * v[353] - v[11] * v[354]
		- v[14] * v[355])*v[74];
	invM(2,4) = v[199] * v[74];
	invM(3,0) = v[74] * (v[22] * v[356] + v[10] * v[358] + v[359] * v[6] - v[341] * v[7] + v[357] * v[8]);
	invM(3,1) = (v[2] * v[341] - v[3] * v[357] - v[1] * v[359] + v[22] * v[360] - v[358] * v[5])*v[74];
	invM(3,2) = (-(v[17] * v[347]) + v[22] * v[361] + v[18] * v[362] + v[20] * v[363] + v[16] * v[364])*v[74];
	invM(3,3) = (v[12] * v[347] - v[13] * v[362] - v[15] * v[363] - v[11] * v[364] + v[22] * v[365])*v[74];
	invM(3,4) = -(v[225] * v[74]);
	invM(4,0) = v[74] * (v[21] * (-v[28] + v[38]) + v[367] * v[6] + v[368] * v[7] - v[348] * v[8] + v[366] * v[9]);
	invM(4,1) = (v[3] * v[348] - v[1] * v[367] - v[2] * v[368] + v[21] * (v[29] - v[39]) - v[366] * v[4])*v[74];
	invM(4,2) = (-(v[18] * v[353]) + v[19] * v[369] + v[21] * (-v[36] + v[37]) + v[16] * v[370] + v[17] * v[371])*v[74];
	invM(4,3) = (v[13] * v[353] - v[14] * v[369] - v[11] * v[370] - v[12] * v[371] + v[21] * (v[46] - v[47]))*v[74];
	invM(4,4) = v[243] * v[74];

	return invM;
}


//Inverte uma matriz de 0x0 a 6x6 escolhendo automaticamente a função correta
Matrix invert(Matrix &matrix)
{
	if (matrix.getLines() == 2 && matrix.getColumns() == 2)
		return invert2x2(matrix);
	else
	{
		if (matrix.getLines() == 3 && matrix.getColumns() == 3)
			return invert3x3(matrix);
		else
		{
			if (matrix.getLines() == 4 && matrix.getColumns() == 4)
				return invert4x4(matrix);
			else
			{
				if (matrix.getLines() == 5 && matrix.getColumns() == 5)
					return invert5x5(matrix);
				else
				{
					if (matrix.getLines() == 6 && matrix.getColumns() == 6)
						return invert6x6(matrix);
					else
					{
						if (matrix.getLines() == 1 && matrix.getColumns() == 1)
						{
							Matrix a(1, 1);
							a(0, 0) = 1.0 / matrix(0, 0);
							return a;
						}
						else
						{
							if (matrix.getLines() == 0 && matrix.getColumns() == 0)
							{
								Matrix a(0, 0);
								return a;
							}
							else
							{
								printf("Uncompatible dimensions. Function 'invert'.\n");
								return 0;
							}
						}
					}
				}
			}
		}
	}
}

//Inverte uma matriz 6x6
Matrix invert6x6(Matrix &matrix)
{
	//Verificação da possibilidade de inversão
	if (matrix.getColumns() != 6 || matrix.getLines() != 6)
	{
		printf("Uncompatible dimensions for input matrix.\n");
		return 0;
	}
	else
	{
		double v[1000];
		Matrix Minv(6, 6);
		
		
		double** pM;
		pM = new double*[6];
		for (int i = 0; i < 6; i++)
			pM[i] = new double[6];
		double** pMinv;
		pMinv = new double*[6];
		for (int i = 0; i < 6; i++)
			pMinv[i] = new double[6];

		matrix.MatrixToPtr(pM, 6);

#pragma region Acegen
		v[976] = -(pM[0][0] * pM[1][2] * pM[3][1]);
		v[975] = pM[0][2] * pM[1][0] * pM[3][1];
		v[974] = pM[0][1] * pM[1][2] * pM[3][0];
		v[973] = -(pM[0][2] * pM[1][1] * pM[3][0]);
		v[966] = -(pM[2][2] * pM[3][1]);
		v[965] = pM[2][0] * pM[3][1];
		v[964] = pM[2][2] * pM[3][0];
		v[963] = -(pM[2][1] * pM[3][0]);
		v[959] = pM[1][0] * v[966];
		v[958] = pM[1][2] * v[965];
		v[957] = pM[1][1] * v[964];
		v[956] = pM[1][2] * v[963];
		v[952] = pM[2][0] * pM[4][1];
		v[951] = -(pM[0][1] * pM[2][5] * pM[4][0]);
		v[945] = pM[1][1] * pM[4][0];
		v[944] = -(pM[0][5] * v[945]);
		v[938] = pM[2][4] * pM[4][0];
		v[937] = pM[2][1] * pM[4][0];
		v[933] = -(pM[0][5] * v[937]);
		v[928] = pM[3][4] * v[952];
		v[927] = pM[1][1] * v[938];
		v[926] = pM[3][4] * pM[4][0];
		v[925] = pM[2][1] * v[926];
		v[921] = pM[0][0] * pM[1][5];
		v[918] = -(pM[2][4] * v[921]);
		v[914] = pM[1][3] * pM[4][5];
		v[911] = pM[1][5] * pM[3][4];
		v[904] = pM[0][3] * pM[2][0];
		v[903] = pM[3][3] * pM[4][5];
		v[902] = pM[0][3] * pM[4][5];
		v[901] = pM[0][0] * pM[4][4];
		v[900] = pM[4][4] * v[904];
		v[899] = pM[3][3] * v[901];
		v[896] = pM[3][5] * pM[4][0];
		v[895] = -(pM[0][4] * pM[2][5]);
		v[884] = pM[1][3] * pM[2][0];
		v[883] = pM[1][0] * v[903];
		v[882] = pM[3][0] * v[914];
		v[887] = -v[882] + v[883];
		v[881] = pM[2][3] * pM[4][5];
		v[880] = pM[1][0] * pM[4][4];
		v[879] = pM[4][4] * v[884];
		v[878] = pM[3][3] * v[880];
		v[877] = pM[3][0] * pM[4][4];
		v[876] = pM[2][3] * v[877];
		v[874] = pM[2][0] * v[911];
		v[873] = -(pM[2][4] * v[896]);
		v[872] = pM[1][4] * pM[4][0];
		v[871] = pM[1][5] * pM[4][0];
		v[870] = -(pM[2][5] * v[872]);
		v[869] = pM[2][4] * v[871];
		v[889] = v[869] + v[870];
		v[859] = pM[0][5] * pM[2][4] * pM[4][2];
		v[850] = pM[0][5] * pM[3][4];
		v[849] = pM[0][2] * pM[4][5];
		v[848] = -(pM[3][4] * v[902]);
		v[847] = pM[1][4] * v[849];
		v[836] = pM[2][1] * v[850];
		v[835] = pM[3][5] * pM[4][1];
		v[834] = pM[4][1] * v[895];
		v[830] = pM[0][4] * pM[4][5];
		v[828] = pM[0][4] * pM[4][3];
		v[827] = -(pM[2][3] * pM[4][2]);
		v[826] = pM[3][3] * pM[4][2];
		v[819] = pM[1][5] * pM[2][1];
		v[815] = pM[3][4] * v[819];
		v[814] = -(pM[2][4] * v[835]);
		v[813] = pM[1][4] * pM[4][1];
		v[812] = pM[1][5] * pM[4][1];
		v[811] = -(pM[2][5] * v[813]);
		v[810] = pM[2][4] * v[812];
		v[823] = v[810] + v[811];
		v[806] = pM[1][2] * pM[4][5];
		v[805] = pM[2][2] * pM[4][5];
		v[804] = pM[2][4] * v[806];
		v[803] = pM[1][4] * v[805];
		v[818] = v[803] - v[804];
		v[802] = pM[3][2] * pM[4][5];
		v[801] = pM[1][4] * v[802];
		v[799] = pM[1][4] * pM[4][3];
		v[798] = pM[3][4] * pM[4][3];
		v[797] = pM[2][2] * v[798];
		v[796] = pM[3][2] * v[799];
		v[795] = pM[2][4] * pM[4][3];
		v[794] = pM[1][3] * pM[3][4];
		v[793] = pM[1][4] * v[826];
		v[765] = pM[2][1] * pM[4][5];
		v[764] = -(pM[3][4] * v[765]);
		v[763] = pM[3][1] * v[881];
		v[762] = pM[0][1] * pM[4][5];
		v[761] = pM[1][1] * v[830];
		v[760] = -(pM[1][4] * v[762]);
		v[759] = pM[0][1] * pM[1][5];
		v[758] = pM[2][5] * pM[4][4];
		v[757] = pM[0][5] * pM[4][4];
		v[756] = -(pM[0][1] * v[758]);
		v[755] = pM[1][1] * v[757];
		v[754] = -(pM[0][5] * pM[1][4]);
		v[753] = pM[0][5] * v[795];
		v[752] = pM[0][4] * pM[1][5];
		v[751] = pM[2][0] * v[828];
		v[750] = -(pM[4][2] * v[759]);
		v[749] = pM[2][5] * pM[4][2];
		v[748] = pM[0][4] * v[749];
		v[747] = pM[3][0] * v[827];
		v[746] = pM[2][3] * pM[3][5];
		v[745] = pM[0][3] * pM[4][2];
		v[744] = pM[1][3] * pM[4][2];
		v[743] = pM[2][1] * pM[4][2];
		v[742] = pM[3][1] * pM[4][2];
		v[741] = pM[2][3] * v[742];
		v[740] = pM[2][0] * pM[3][3];
		v[739] = -(pM[3][3] * v[743]);
		v[738] = pM[1][1] * v[745];
		v[737] = -(pM[0][1] * v[744]);
		v[736] = -(pM[3][0] * pM[4][1]);
		v[735] = pM[4][1] * v[740];
		v[734] = pM[4][1] * v[746];
		v[733] = -(pM[1][0] * pM[4][1]);
		v[732] = pM[0][0] * pM[4][1];
		v[731] = pM[0][3] * pM[1][5];
		v[730] = pM[2][2] * v[736];
		v[729] = pM[3][2] * pM[4][1];
		v[728] = -(pM[2][3] * v[729]);
		v[727] = pM[0][2] * v[733];
		v[726] = pM[1][2] * v[732];
		v[725] = pM[2][1] * pM[3][5];
		v[724] = pM[0][2] * v[725];
		v[723] = pM[1][1] * pM[3][5];
		v[722] = pM[0][1] * pM[3][5];
		v[721] = -(pM[0][2] * v[723]);
		v[720] = pM[1][2] * v[722];
		v[980] = -v[720] - v[721];
		v[719] = pM[2][0] * pM[3][5];
		v[718] = pM[0][2] * v[719];
		v[717] = pM[1][0] * pM[3][5];
		v[716] = pM[0][0] * pM[3][5];
		v[715] = -(pM[0][2] * v[717]);
		v[714] = pM[1][2] * v[716];
		v[713] = pM[0][5] * v[794];
		v[712] = -(pM[3][4] * v[731]);
		v[711] = pM[2][3] * pM[3][4];
		v[710] = pM[1][5] * v[711];
		v[709] = pM[2][5] * pM[3][4];
		v[708] = -(pM[1][3] * v[709]);
		v[707] = pM[3][3] * v[754];
		v[706] = pM[3][3] * v[752];
		v[917] = -v[706] - v[707] - v[713];
		v[705] = pM[2][4] * pM[3][3];
		v[704] = -(pM[1][5] * v[705]);
		v[703] = pM[2][5] * pM[3][3];
		v[702] = pM[1][4] * v[703];
		v[701] = pM[2][2] * pM[3][3];
		v[700] = -(pM[1][2] * pM[3][3]);
		v[699] = pM[1][1] * v[701];
		v[698] = pM[2][1] * v[700];
		v[697] = pM[0][1] * pM[3][3];
		v[696] = pM[0][2] * pM[3][3];
		v[695] = pM[1][2] * v[697];
		v[694] = -(pM[1][1] * v[696]);
		v[693] = pM[2][1] * pM[3][2];
		v[692] = -(pM[1][1] * pM[3][2]);
		v[691] = pM[0][5] * pM[3][2];
		v[690] = pM[1][5] * pM[3][2];
		v[689] = pM[1][0] * v[691];
		v[688] = -(pM[0][0] * v[690]);
		v[687] = pM[2][5] * pM[3][2];
		v[686] = pM[0][0] * v[687];
		v[685] = pM[0][4] * pM[1][3];
		v[684] = pM[0][3] * pM[1][4];
		v[683] = pM[0][4] * pM[2][3];
		v[682] = pM[1][4] * pM[2][3];
		v[681] = pM[0][3] * pM[2][4];
		v[680] = pM[1][3] * pM[2][4];
		v[679] = pM[2][3] * v[692];
		v[678] = pM[1][3] * v[693];
		v[677] = pM[0][1] * pM[3][2];
		v[676] = pM[0][3] * pM[3][2];
		v[675] = -(pM[1][3] * v[677]);
		v[674] = pM[1][1] * v[676];
		v[673] = pM[0][5] * pM[3][1];
		v[672] = pM[1][5] * pM[3][1];
		v[671] = pM[1][0] * v[673];
		v[670] = -(pM[0][0] * v[672]);
		v[669] = pM[2][5] * pM[3][1];
		v[668] = pM[0][0] * v[669];
		v[667] = pM[0][4] * pM[3][1];
		v[666] = pM[1][4] * pM[3][1];
		v[665] = pM[1][0] * v[667];
		v[664] = -(pM[0][0] * v[666]);
		v[663] = pM[2][4] * pM[3][1];
		v[662] = pM[0][0] * v[663];
		v[661] = pM[0][5] * pM[3][0];
		v[660] = pM[1][5] * pM[3][0];
		v[659] = -(pM[1][1] * v[661]);
		v[658] = pM[0][1] * v[660];
		v[657] = pM[2][5] * pM[3][0];
		v[656] = -(pM[0][1] * v[657]);
		v[655] = pM[0][4] * pM[3][0];
		v[654] = pM[1][4] * pM[3][0];
		v[653] = -(pM[1][1] * v[655]);
		v[652] = pM[0][1] * v[654];
		v[651] = pM[2][4] * pM[3][0];
		v[650] = -(pM[0][1] * v[651]);
		v[649] = pM[2][5] * v[685];
		v[648] = -(pM[2][5] * v[684]);
		v[647] = pM[0][2] * pM[2][5];
		v[646] = -(pM[1][2] * pM[2][5]);
		v[645] = pM[1][0] * v[647];
		v[644] = pM[0][0] * v[646];
		v[643] = pM[0][1] * pM[1][0];
		v[642] = -(pM[0][0] * pM[1][1]);
		v[641] = -(pM[0][5] * v[680]);
		v[640] = pM[1][5] * v[681];
		v[639] = pM[0][3] * pM[1][2];
		v[638] = pM[0][2] * pM[1][3];
		v[637] = pM[2][4] * v[643];
		v[636] = pM[2][4] * v[642];
		v[635] = pM[0][5] * v[682];
		v[634] = -(pM[1][5] * v[683]);
		v[862] = v[634] + v[635] + v[641];
		v[767] = v[640] + v[648] + v[649] + v[862];
		v[633] = pM[1][2] * pM[2][3];
		v[632] = pM[0][2] * pM[2][3];
		v[631] = pM[0][1] * v[633];
		v[630] = pM[1][1] * v[632];
		v[629] = -(pM[0][5] * pM[2][2]);
		v[628] = pM[1][5] * pM[2][2];
		v[627] = pM[1][0] * v[629];
		v[626] = pM[0][0] * v[628];
		v[625] = pM[1][3] * pM[2][2];
		v[624] = pM[0][3] * pM[2][2];
		v[623] = pM[0][1] * v[625];
		v[622] = pM[1][1] * v[624];
		v[621] = pM[0][5] * pM[1][2];
		v[620] = -(pM[0][2] * pM[1][5]);
		v[619] = -(pM[1][0] * pM[2][1]);
		v[618] = pM[0][0] * pM[2][1];
		v[617] = pM[0][4] * v[619];
		v[616] = pM[1][4] * v[618];
		v[615] = pM[2][1] * v[639];
		v[614] = pM[2][1] * v[638];
		v[844] = -v[614] + v[615] - v[622] + v[623] + v[630] - v[631];
		v[613] = pM[2][0] * v[621];
		v[612] = pM[2][0] * v[620];
		v[775] = v[612] + v[613] + v[626] + v[627] + v[644] + v[645];
		v[611] = pM[1][1] * pM[2][0];
		v[610] = -(pM[0][1] * pM[2][0]);
		v[609] = pM[0][4] * v[611];
		v[608] = pM[1][4] * v[610];
		v[790] = v[608] + v[609] + v[616] + v[617] + v[636] + v[637];
		v[156] = pM[1][5] * v[610];
		v[155] = pM[0][5] * v[611];
		v[158] = pM[1][5] * v[618];
		v[157] = pM[0][5] * v[619];
		v[77] = pM[2][1] * v[620];
		v[76] = pM[2][1] * v[621];
		v[118] = -(pM[1][4] * v[624]);
		v[117] = pM[0][4] * v[625];
		v[79] = pM[0][1] * v[628];
		v[78] = pM[1][1] * v[629];
		v[120] = pM[1][4] * v[632];
		v[119] = -(pM[0][4] * v[633]);
		v[122] = -(pM[2][4] * v[638]);
		v[121] = pM[2][4] * v[639];
		v[778] = v[117] + v[118] + v[119] + v[120] + v[121] + v[122];
		v[160] = pM[2][5] * v[642];
		v[159] = pM[2][5] * v[643];
		v[782] = v[155] + v[156] + v[157] + v[158] + v[159] + v[160];
		v[81] = pM[0][1] * v[646];
		v[80] = pM[1][1] * v[647];
		v[770] = v[76] + v[77] + v[78] + v[79] + v[80] + v[81];
		v[174] = pM[1][1] * v[651];
		v[171] = -(pM[2][1] * v[654]);
		v[170] = pM[2][1] * v[655];
		v[147] = pM[1][1] * v[657];
		v[144] = -(pM[2][1] * v[660]);
		v[143] = pM[2][1] * v[661];
		v[97] = -(pM[3][0] * v[647]);
		v[96] = -(pM[3][0] * v[646]);
		v[95] = -(pM[3][0] * v[620]);
		v[979] = v[714] + v[715] + v[95];
		v[94] = -(pM[3][0] * v[621]);
		v[776] = v[688] + v[689] + v[94] + v[979];
		v[93] = -(pM[3][0] * v[628]);
		v[92] = -(pM[3][0] * v[629]);
		v[313] = pM[3][1] * v[638];
		v[311] = -(pM[3][1] * v[639]);
		v[846] = v[311] + v[313] + v[674] + v[675] + v[694] + v[695];
		v[268] = pM[3][1] * v[632];
		v[266] = pM[3][1] * v[624];
		v[219] = pM[3][1] * v[633];
		v[217] = -(pM[3][1] * v[625]);
		v[800] = v[217] + v[219] + v[678] + v[679] + v[698] + v[699];
		v[180] = -(pM[1][0] * v[663]);
		v[177] = pM[2][0] * v[666];
		v[176] = -(pM[2][0] * v[667]);
		v[153] = -(pM[1][0] * v[669]);
		v[150] = pM[2][0] * v[672];
		v[149] = -(pM[2][0] * v[673]);
		v[69] = -(pM[3][1] * v[647]);
		v[68] = -(pM[3][1] * v[646]);
		v[67] = -(pM[3][1] * v[620]);
		v[66] = -(pM[3][1] * v[621]);
		v[65] = -(pM[3][1] * v[628]);
		v[64] = -(pM[3][1] * v[629]);
		v[260] = pM[2][1] * v[676];
		v[258] = pM[2][3] * v[677];
		v[128] = pM[3][2] * v[680];
		v[127] = -(pM[3][2] * v[681]);
		v[126] = -(pM[3][2] * v[682]);
		v[125] = pM[3][2] * v[683];
		v[124] = pM[3][2] * v[684];
		v[123] = -(pM[3][2] * v[685]);
		v[102] = -(pM[1][0] * v[687]);
		v[99] = pM[2][0] * v[690];
		v[98] = -(pM[2][0] * v[691]);
		v[75] = pM[2][5] * v[677];
		v[74] = pM[2][5] * v[692];
		v[73] = -(pM[1][5] * v[677]);
		v[72] = pM[1][1] * v[691];
		v[787] = v[72] + v[720] + v[721] + v[73];
		v[772] = v[66] + v[67] + v[787];
		v[71] = pM[1][5] * v[693];
		v[70] = -(pM[2][1] * v[691]);
		v[264] = pM[2][1] * v[696];
		v[262] = pM[2][2] * v[697];
		v[829] = v[258] - v[260] - v[262] + v[264] + v[266] - v[268];
		v[134] = pM[2][4] * v[700];
		v[133] = pM[2][4] * v[696];
		v[132] = pM[1][4] * v[701];
		v[131] = -(pM[0][4] * v[701]);
		v[130] = -(pM[1][4] * v[696]);
		v[129] = -(pM[0][4] * v[700]);
		v[48] = -(pM[0][4] * v[703]);
		v[46] = pM[0][5] * v[705];
		v[193] = -(pM[3][4] * v[642]);
		v[192] = -(pM[3][4] * v[643]);
		v[948] = v[192] + v[193] + v[664] + v[665];
		v[792] = v[652] + v[653] + v[948];
		v[191] = -(pM[3][4] * v[618]);
		v[190] = -(pM[3][4] * v[619]);
		v[189] = -(pM[3][4] * v[610]);
		v[791] = v[170] + v[176] + v[189] + v[191] + v[650] + v[662];
		v[188] = -(pM[3][4] * v[611]);
		v[929] = v[171] + v[174] + v[177] + v[188];
		v[789] = v[180] + v[190] + v[929];
		v[140] = pM[3][4] * v[633];
		v[139] = -(pM[3][4] * v[632]);
		v[138] = -(pM[3][4] * v[625]);
		v[781] = v[126] + v[128] + v[132] + v[134] + v[138] + v[140];
		v[137] = pM[3][4] * v[624];
		v[780] = v[125] + v[127] + v[131] + v[133] + v[137] + v[139];
		v[136] = pM[3][4] * v[638];
		v[135] = -(pM[3][4] * v[639]);
		v[779] = v[123] + v[124] + v[129] + v[130] + v[135] + v[136];
		v[54] = pM[0][3] * v[709];
		v[52] = -(pM[0][5] * v[711]);
		v[841] = v[46] + v[52];
		v[166] = -(pM[3][5] * v[642]);
		v[165] = -(pM[3][5] * v[643]);
		v[978] = -v[165] - v[166] - v[670];
		v[785] = v[165] + v[166] + v[658] + v[659] + v[670] + v[671];
		v[164] = -(pM[3][5] * v[618]);
		v[163] = -(pM[3][5] * v[619]);
		v[962] = -v[144] - v[150] - v[153] - v[163];
		v[162] = -(pM[3][5] * v[610]);
		v[941] = v[162] + v[164] + v[668];
		v[784] = v[143] + v[149] + v[656] + v[941];
		v[161] = -(pM[3][5] * v[611]);
		v[783] = v[144] + v[147] + v[150] + v[153] + v[161] + v[163];
		v[113] = -(pM[2][2] * v[716]);
		v[777] = v[113] + v[686] + v[718] + v[92] + v[97] + v[98];
		v[112] = pM[2][2] * v[717];
		v[961] = v[102] + v[112] + v[93] + v[99];
		v[110] = -(pM[1][2] * v[719]);
		v[774] = v[110] + v[96] + v[961];
		v[85] = -(pM[2][2] * v[722]);
		v[972] = v[724] + v[85];
		v[788] = v[70] + v[75] + v[972];
		v[773] = v[64] + v[69] + v[788];
		v[84] = pM[2][2] * v[723];
		v[82] = -(pM[1][2] * v[725]);
		v[786] = v[71] + v[74] + v[82] + v[84];
		v[771] = v[65] + v[68] + v[786];
		v[61] = pM[3][5] * v[680];
		v[60] = -(pM[3][5] * v[681]);
		v[59] = -(pM[3][5] * v[682]);
		v[768] = v[59] + v[61] + v[702] + v[704] + v[708] + v[710];
		v[58] = pM[3][5] * v[683];
		v[766] = v[48] + v[54] + v[58] + v[60] + v[841];
		v[57] = pM[3][5] * v[684];
		v[56] = -(pM[3][5] * v[685]);
		v[769] = v[56] + v[57] + v[706] + v[707] + v[712] + v[713];
		v[569] = pM[4][1] * v[638];
		v[568] = -(pM[4][1] * v[639]);
		v[982] = v[568] + v[569] + v[737] + v[738];
		v[552] = pM[4][1] * v[701];
		v[968] = v[552] + v[728] + v[739] + v[741];
		v[550] = pM[2][0] * v[729];
		v[492] = pM[4][1] * v[731];
		v[491] = pM[1][3] * v[732];
		v[490] = pM[0][3] * v[733];
		v[953] = v[490] + v[491];
		v[465] = -(pM[4][1] * v[703]);
		v[935] = v[465] + v[734];
		v[461] = -(pM[4][1] * v[719]);
		v[458] = pM[4][1] * v[657];
		v[934] = v[458] + v[461];
		v[457] = pM[2][3] * v[736];
		v[936] = v[457] + v[735];
		v[575] = pM[4][2] * v[642];
		v[574] = pM[4][2] * v[643];
		v[981] = v[574] + v[575] + v[726] + v[727];
		v[556] = pM[4][2] * v[740];
		v[969] = v[556] + v[747];
		v[554] = -(pM[2][0] * v[742]);
		v[553] = pM[3][0] * v[743];
		v[967] = v[550] + v[553] + v[554] + v[730];
		v[428] = pM[0][0] * v[744];
		v[427] = -(pM[1][0] * v[745]);
		v[922] = v[427] + v[428];
		v[410] = pM[4][2] * v[746];
		v[409] = -(pM[4][2] * v[703]);
		v[905] = v[409] + v[410];
		v[398] = -(pM[4][2] * v[719]);
		v[373] = pM[1][4] * v[749];
		v[875] = v[373] + v[804];
		v[919] = -v[803] + v[875];
		v[413] = pM[4][3] * v[752];
		v[358] = pM[4][3] * v[754];
		v[923] = v[358] + v[413];
		v[463] = pM[2][1] * v[757];
		v[943] = v[463] + v[756] - v[834];
		v[447] = pM[1][1] * v[758];
		v[360] = pM[4][4] * v[759];
		v[863] = v[360] + v[760] + v[761];
		v[949] = -v[755] + v[863];
		v[359] = -(pM[4][4] * v[731]);
		v[860] = v[358] + v[359];
		v[326] = -(pM[3][3] * v[757]);
		v[325] = pM[4][4] * v[673];
		v[353] = pM[1][3] * v[762];
		v[277] = pM[4][5] * v[711];
		v[274] = -(pM[4][5] * v[705]);
		v[833] = v[274] + v[277];
		v[273] = pM[3][3] * v[765];
		v[840] = v[273] - v[763];
		v[271] = pM[4][5] * v[663];
		v[942] = -v[271] - v[764];
		v[37] = pM[1][2] * v[766] + pM[3][2] * v[767] + pM[0][2] * v[768] + pM[2][2] * v[769];
		v[62] = pM[1][1] * v[766] + pM[3][1] * v[767] + pM[0][1] * v[768] + pM[2][1] * v[769];
		v[63] = pM[3][4] * v[770] + pM[0][4] * v[771] + pM[2][4] * v[772] + pM[1][4] * v[773];
		v[88] = pM[3][3] * v[770] + pM[0][3] * v[771] + pM[2][3] * v[772] + pM[1][3] * v[773];
		v[89] = pM[3][1] * v[778] + pM[2][1] * v[779] + pM[1][1] * v[780] + pM[0][1] * v[781];
		v[363] = pM[4][1] * v[37] - pM[4][2] * v[62] + pM[4][3] * v[63] - pM[4][4] * v[88] + pM[4][5] * v[89];
		v[90] = pM[1][0] * v[766] + pM[3][0] * v[767] + pM[0][0] * v[768] + pM[2][0] * v[769];
		v[91] = pM[0][4] * v[774] + pM[3][4] * v[775] + pM[2][4] * v[776] + pM[1][4] * v[777];
		v[116] = pM[0][3] * v[774] + pM[3][3] * v[775] + pM[2][3] * v[776] + pM[1][3] * v[777];
		v[141] = pM[3][0] * v[778] + pM[2][0] * v[779] + pM[1][0] * v[780] + pM[0][0] * v[781];
		v[431] = -(pM[4][4] * v[116]) + pM[4][5] * v[141] + pM[4][0] * v[37] - pM[4][2] * v[90] + pM[4][3] * v[91];
		v[142] = pM[3][4] * v[782] + pM[0][4] * v[783] + pM[1][4] * v[784] + pM[2][4] * v[785];
		v[167] = pM[3][3] * v[782] + pM[0][3] * v[783] + pM[1][3] * v[784] + pM[2][3] * v[785];
		v[168] = pM[0][3] * v[789] + pM[3][3] * v[790] + pM[1][3] * v[791] + pM[2][3] * v[792];
		v[495] = pM[4][3] * v[142] - pM[4][4] * v[167] + pM[4][5] * v[168] + pM[4][0] * v[62] - pM[4][1] * v[90];
		v[169] = -(pM[3][1] * v[775]) + pM[0][0] * v[786] + pM[2][0] * v[787] + pM[1][0] * v[788] + pM[2][1] * (-v[94] - v[95])
			+ pM[0][1] * (-v[93] - v[96]) + pM[1][1] * (-v[92] - v[97]);
		v[194] = pM[0][2] * v[789] + pM[3][2] * v[790] + pM[1][2] * v[791] + pM[2][2] * v[792];
		v[538] = pM[4][2] * v[142] - pM[4][4] * v[169] + pM[4][5] * v[194] + pM[4][0] * v[63] - pM[4][1] * v[91];
		v[195] = pM[0][0] * v[800] + pM[1][0] * v[829] + pM[3][0] * v[844] + pM[2][0] * v[846];
		v[605] = -(pM[4][1] * v[141]) + pM[4][2] * v[168] - pM[4][3] * v[194] + pM[4][4] * v[195] + pM[4][0] * v[89];
		v[578] = -(pM[4][1] * v[116]) + pM[4][2] * v[167] - pM[4][3] * v[169] + pM[4][5] * v[195] + pM[4][0] * v[88];
		v[196] = 1e0 / (-(pM[5][0] * v[363]) + pM[5][1] * v[431] - pM[5][2] * v[495] + pM[5][3] * v[538] - pM[5][4] * v[578]
			+ pM[5][5] * v[605]);
		v[197] = pM[4][2] * v[705];
		v[199] = pM[4][2] * v[711];
		v[200] = pM[4][2] * v[794];
		v[201] = pM[4][2] * v[682];
		v[202] = pM[4][2] * v[680];
		v[203] = pM[3][2] * v[795];
		v[809] = -v[197] + v[199] + v[203] - v[797];
		v[206] = pM[1][2] * v[798];
		v[808] = -v[200] + v[206] + v[793] - v[796];
		v[207] = pM[2][2] * v[799];
		v[208] = pM[1][2] * v[795];
		v[807] = -v[201] + v[202] + v[207] - v[208];
		v[580] = -(pM[4][1] * v[781]) + pM[4][4] * v[800] + pM[3][1] * v[807] + pM[2][1] * v[808] + pM[1][1] * v[809];
		v[210] = pM[2][3] * v[690];
		v[212] = -(pM[1][3] * v[687]);
		v[214] = -(pM[3][3] * v[628]);
		v[216] = -(pM[3][3] * v[646]);
		v[218] = pM[3][5] * v[625];
		v[220] = -(pM[3][5] * v[633]);
		v[820] = v[210] + v[212] + v[214] + v[216] + v[218] + v[220];
		v[222] = pM[2][4] * v[802];
		v[225] = pM[3][4] * v[805];
		v[816] = v[222] - v[225];
		v[226] = pM[3][4] * v[806];
		v[817] = v[226] - v[801];
		v[365] = -(pM[3][5] * v[807]) - pM[2][5] * v[808] - pM[1][5] * v[809] + pM[1][3] * v[816] + pM[2][3] * v[817]
			+ pM[3][3] * v[818] + pM[4][4] * v[820];
		v[229] = -(pM[3][4] * v[812]);
		v[230] = pM[4][1] * v[709];
		v[824] = v[230] + v[814];
		v[838] = v[271] + v[764] + v[824];
		v[231] = pM[3][5] * v[813];
		v[821] = v[229] + v[231];
		v[233] = -(pM[1][5] * v[663]);
		v[234] = pM[2][5] * v[666];
		v[236] = -(pM[1][1] * v[709]);
		v[237] = -(pM[1][4] * v[725]);
		v[238] = pM[2][4] * v[723];
		v[825] = v[233] + v[234] + v[236] + v[237] + v[238] + v[815];
		v[497] = -(pM[4][4] * v[771]) - pM[1][1] * v[816] - pM[2][1] * v[817] - pM[3][1] * v[818] + pM[2][2] * v[821]
			+ pM[3][2] * v[823] + pM[1][2] * v[824] + pM[4][2] * v[825];
		v[239] = pM[2][3] * v[672];
		v[240] = pM[1][3] * v[669];
		v[241] = pM[3][3] * v[819];
		v[242] = pM[1][1] * v[703];
		v[243] = pM[1][3] * v[725];
		v[244] = pM[2][3] * v[723];
		v[822] = -v[239] + v[240] + v[241] - v[242] - v[243] + v[244];
		v[540] = -(pM[4][3] * v[771]) + pM[4][5] * v[800] + pM[4][1] * v[820] + pM[4][2] * v[822];
		v[433] = pM[2][3] * v[821] - pM[4][4] * v[822] + pM[3][3] * v[823] + pM[4][3] * v[825] + pM[1][1] * v[833]
			+ pM[1][3] * v[838] + pM[1][4] * v[840];
		v[246] = pM[0][5] * v[826];
		v[247] = pM[0][4] * v[826];
		v[248] = pM[0][5] * v[827];
		v[249] = pM[3][4] * v[745];
		v[250] = pM[4][2] * v[683];
		v[251] = pM[4][2] * v[681];
		v[252] = pM[2][4] * v[691];
		v[253] = pM[3][2] * v[828];
		v[254] = -(pM[3][4] * v[629]);
		v[837] = -v[252] + v[254];
		v[255] = pM[0][2] * v[798];
		v[831] = v[247] - v[249] - v[253] + v[255];
		v[256] = pM[2][2] * v[828];
		v[257] = pM[0][2] * v[795];
		v[832] = -v[250] + v[251] + v[256] - v[257];
		v[586] = pM[4][1] * v[780] + pM[0][1] * v[809] - pM[4][4] * v[829] + pM[2][1] * v[831] + pM[3][1] * v[832];
		v[259] = pM[2][3] * v[691];
		v[261] = -(pM[2][5] * v[676]);
		v[263] = pM[3][3] * v[629];
		v[265] = pM[3][3] * v[647];
		v[267] = pM[3][5] * v[624];
		v[269] = -(pM[3][5] * v[632]);
		v[839] = v[259] + v[261] + v[263] + v[265] + v[267] + v[269];
		v[270] = pM[3][2] * v[830];
		v[272] = pM[0][4] * v[805];
		v[391] = pM[2][4] * v[246] + pM[3][4] * v[248] - pM[2][3] * v[270] + pM[3][3] * v[272] + pM[0][3] * v[816]
			- pM[2][5] * v[831] - pM[3][5] * v[832] + pM[0][2] * v[833] + pM[4][3] * v[837] + pM[4][4] * v[839];
		v[279] = pM[0][4] * v[835];
		v[852] = v[279] + v[325];
		v[280] = -(pM[0][5] * v[663]);
		v[281] = pM[2][5] * v[667];
		v[283] = -(pM[0][1] * v[709]);
		v[284] = -(pM[0][4] * v[725]);
		v[285] = pM[2][4] * v[722];
		v[843] = v[280] + v[281] + v[283] + v[284] + v[285] + v[836];
		v[507] = pM[2][1] * v[270] - pM[3][1] * v[272] + pM[2][2] * v[279] + pM[4][4] * v[773] - pM[0][1] * v[816]
			+ pM[3][2] * v[834] - pM[4][1] * v[837] + pM[0][2] * v[838] + pM[4][2] * v[843];
		v[286] = pM[0][3] * v[669];
		v[287] = pM[2][5] * v[697];
		v[288] = pM[0][3] * v[725];
		v[289] = pM[2][3] * v[722];
		v[842] = v[286] - v[287] - v[288] + v[289];
		v[545] = pM[2][1] * v[246] + pM[3][1] * v[248] + pM[4][5] * (-v[258] + v[260] + v[262] - v[266]) + pM[4][3] * v[773]
			+ pM[4][1] * v[839] - pM[0][2] * v[840] + pM[4][2] * v[842];
		v[449] = pM[2][1] * v[326] + pM[0][1] * v[833] + pM[3][3] * v[834] + pM[0][3] * v[838] + pM[0][4] * v[840]
			+ pM[4][1] * v[841] - pM[4][4] * v[842] + pM[4][3] * v[843] + pM[2][3] * v[852];
		v[291] = -(pM[4][2] * v[754]);
		v[292] = -(pM[1][5] * v[826]);
		v[293] = -(pM[0][5] * v[744]);
		v[294] = pM[4][2] * v[731];
		v[854] = v[293] + v[294];
		v[295] = pM[4][2] * v[685];
		v[296] = pM[4][2] * v[684];
		v[297] = pM[1][4] * v[691];
		v[298] = pM[0][4] * v[690];
		v[299] = pM[3][4] * v[621];
		v[851] = -v[297] + v[298] + v[299];
		v[300] = -(pM[1][5] * v[798]);
		v[301] = pM[1][2] * v[828];
		v[302] = pM[0][2] * v[799];
		v[845] = -v[295] + v[296] + v[301] - v[302];
		v[598] = pM[4][1] * v[778] + pM[0][1] * v[807] - pM[1][1] * v[832] - pM[4][4] * v[844] + pM[2][1] * v[845];
		v[592] = -(pM[4][1] * v[779]) - pM[0][1] * v[808] + pM[1][1] * v[831] + pM[3][1] * v[845] + pM[4][4] * v[846];
		v[304] = pM[1][3] * v[691];
		v[306] = -(pM[1][5] * v[676]);
		v[308] = -(pM[3][3] * v[621]);
		v[310] = -(pM[3][3] * v[620]);
		v[312] = pM[3][5] * v[639];
		v[314] = -(pM[3][5] * v[638]);
		v[855] = v[304] + v[306] + v[308] + v[310] + v[312] + v[314];
		v[315] = -(pM[4][5] * v[685]);
		v[316] = pM[4][5] * v[684];
		v[858] = v[315] + v[316];
		v[861] = v[858] + v[860];
		v[317] = pM[0][4] * v[806];
		v[853] = v[291] + v[317] - v[847];
		v[320] = pM[3][4] * v[849];
		v[897] = v[270] - v[320];
		v[412] = pM[0][4] * v[292] + pM[0][2] * v[300] + pM[1][3] * v[320] - pM[3][5] * v[845] + pM[1][2] * v[848]
			+ pM[4][3] * v[851] + pM[3][3] * v[853] + pM[3][4] * v[854] + pM[4][4] * v[855] + pM[3][2] * v[858];
		v[321] = pM[1][5] * v[667];
		v[322] = pM[1][1] * v[850];
		v[323] = -(pM[0][4] * v[723]);
		v[324] = pM[1][4] * v[722];
		v[856] = v[321] + v[322] + v[323] + v[324];
		v[517] = pM[3][4] * v[750] + pM[0][1] * v[817] - pM[0][2] * v[821] - pM[4][1] * v[851] + pM[1][2] * v[852]
			- pM[3][1] * v[853] + pM[4][2] * v[856] + pM[1][1] * v[897] + pM[4][4] * (-v[67] - v[72] - v[73] + v[980]);
		v[327] = pM[0][3] * v[723];
		v[328] = pM[1][3] * v[722];
		v[857] = -v[327] + v[328];
		v[559] = pM[1][1] * v[246] + pM[0][1] * v[292] - pM[4][3] * v[772] + pM[4][5] * v[846] + pM[3][1] * v[854]
			+ pM[4][1] * v[855] + pM[4][2] * v[857];
		v[468] = pM[0][1] * v[300] + pM[3][4] * v[353] - pM[0][3] * v[821] + pM[1][1] * (v[326] + v[848]) + pM[1][3] * v[852]
			+ pM[4][3] * v[856] - pM[4][4] * v[857] + pM[3][1] * v[861] + pM[3][3] * v[863] + pM[4][1] * v[917];
		v[331] = pM[4][2] * v[752];
		v[868] = -v[291] - v[317] + v[331] + v[847];
		v[333] = pM[2][5] * v[744];
		v[334] = -(pM[2][5] * v[745]);
		v[335] = -(pM[1][4] * v[629]);
		v[336] = pM[0][4] * v[628];
		v[337] = pM[2][4] * v[621];
		v[864] = v[335] - v[336] - v[337];
		v[338] = -(pM[1][5] * v[795]);
		v[339] = -(pM[2][5] * v[828]);
		v[924] = -v[339] - v[753];
		v[340] = pM[2][5] * v[799];
		v[865] = v[338] + v[340];
		v[342] = pM[0][5] * v[625];
		v[344] = -(pM[1][5] * v[624]);
		v[346] = -(pM[2][3] * v[621]);
		v[348] = -(pM[2][3] * v[620]);
		v[350] = pM[2][5] * v[639];
		v[352] = -(pM[2][5] * v[638]);
		v[866] = v[342] + v[344] + v[346] + v[348] + v[350] + v[352];
		v[564] = pM[0][1] * v[333] + pM[1][1] * (-v[248] + v[334]) - pM[2][2] * v[353] + pM[4][5] * (v[614] - v[615] + v[622]
			- v[630] + v[631]) + pM[2][3] * v[750] + pM[4][3] * v[770] + pM[2][1] * v[854] + pM[4][1] * v[866];
		v[356] = -(pM[4][5] * v[681]);
		v[867] = v[339] + v[356];
		v[483] = pM[2][4] * v[353] + pM[0][3] * (v[447] + v[823]) + pM[2][1] * (v[413] + v[861]) + pM[4][1] * v[862]
			+ pM[0][1] * v[865] + pM[1][1] * (v[753] + v[867]) + pM[1][3] * v[943] + pM[2][3] * v[949];
		v[357] = pM[2][4] * v[849];
		v[898] = v[357] - v[859];
		v[920] = -v[272] + v[748] + v[898];
		v[527] = pM[2][4] * v[750] + pM[4][4] * (v[770] - v[79]) + pM[0][2] * v[823] - pM[1][2] * v[834] + pM[1][1] * (-v[357]
			- v[748] + v[859]) + pM[2][2] * v[863] + pM[4][1] * v[864] + pM[2][1] * v[868] + pM[0][1] * v[875];
		v[423] = pM[2][4] * v[294] + pM[0][4] * v[333] + pM[1][4] * v[334] + pM[2][2] * v[858] - pM[4][3] * v[864]
			+ pM[0][2] * v[865] + pM[4][4] * v[866] + pM[1][2] * v[867] - pM[2][3] * v[868] + pM[1][3] * v[898];
		v[368] = -(pM[3][4] * v[871]);
		v[369] = pM[4][0] * v[709];
		v[890] = v[369] + v[873];
		v[370] = pM[3][5] * v[872];
		v[888] = v[368] + v[370];
		v[372] = -(pM[1][5] * v[651]);
		v[375] = -(pM[1][0] * v[709]);
		v[376] = -(pM[1][4] * v[719]);
		v[377] = pM[2][4] * v[717];
		v[891] = v[372] + v[375] + v[376] + v[377] + v[874];
		v[498] = -(pM[4][4] * v[774]) - pM[1][0] * v[816] - pM[2][0] * v[817] + pM[2][2] * v[888] + pM[3][2] * v[889]
			+ pM[1][2] * v[890] + pM[4][2] * v[891] + pM[3][0] * v[919];
		v[379] = -(pM[1][3] * v[877]);
		v[893] = v[379] + v[878];
		v[380] = -(pM[4][4] * v[740]);
		v[892] = v[380] + v[876];
		v[383] = -(pM[2][3] * v[880]);
		v[894] = v[383] + v[879];
		v[581] = -(pM[4][0] * v[781]) + pM[3][0] * v[807] + pM[2][0] * v[808] + pM[1][0] * v[809] + pM[1][2] * v[892]
			+ pM[2][2] * v[893] + pM[3][2] * v[894];
		v[384] = pM[3][0] * v[881];
		v[386] = pM[4][5] * v[740];
		v[885] = v[384] - v[386];
		v[388] = pM[4][5] * v[884];
		v[389] = pM[1][0] * v[881];
		v[886] = v[388] - v[389];
		v[542] = -(pM[4][3] * v[783]) - pM[4][0] * v[822] + pM[1][1] * v[885] + pM[3][1] * v[886] + pM[2][1] * v[887]
			+ pM[1][3] * v[934] + pM[1][0] * v[935] + pM[1][5] * v[936];
		v[541] = -(pM[2][0] * v[292]) + pM[3][0] * v[333] + pM[1][3] * v[398] + pM[1][5] * v[747] - pM[4][3] * v[774]
			+ pM[4][0] * v[820] + pM[1][2] * v[885] + pM[3][2] * v[886] + pM[2][2] * v[887] + pM[1][0] * v[905];
		v[434] = pM[3][0] * v[340] - pM[1][4] * v[885] - pM[3][4] * v[886] - pM[2][4] * v[887] + pM[2][3] * v[888]
			+ pM[3][3] * v[889] + pM[1][3] * v[890] + pM[4][3] * v[891] + pM[1][5] * v[892] + pM[2][5] * v[893] + pM[3][5] * v[894];
		v[392] = pM[4][0] * v[895];
		v[393] = pM[0][4] * v[896];
		v[396] = pM[2][0] * v[850];
		v[397] = -(pM[0][0] * v[709]);
		v[399] = pM[2][4] * v[716];
		v[906] = v[396] + v[397] + v[399];
		v[508] = pM[3][2] * v[392] + pM[2][2] * v[393] + pM[0][4] * v[398] + pM[4][4] * v[777] - pM[0][0] * v[816]
			- pM[4][0] * v[837] + pM[0][2] * v[890] + pM[2][0] * v[897] + pM[4][2] * v[906] + pM[3][0] * v[920];
		v[401] = -(pM[0][3] * v[877]);
		v[909] = v[401] + v[899];
		v[404] = -(pM[2][3] * v[901]);
		v[910] = v[404] + v[900];
		v[912] = -v[751] + v[910];
		v[587] = pM[4][0] * v[780] + pM[0][0] * v[809] + pM[2][0] * v[831] + pM[3][0] * v[832] + pM[0][2] * v[892]
			+ pM[2][2] * v[909] + pM[3][2] * v[910];
		v[405] = pM[3][0] * v[902];
		v[406] = pM[0][0] * v[903];
		v[907] = -v[405] + v[406];
		v[407] = pM[4][5] * v[904];
		v[408] = pM[0][0] * v[881];
		v[908] = v[407] - v[408];
		v[546] = pM[2][0] * v[246] - pM[3][0] * v[334] + pM[0][3] * v[398] + pM[0][5] * v[747] + pM[4][3] * v[777]
			+ pM[4][0] * v[839] + pM[0][2] * v[885] + pM[0][0] * v[905] + pM[2][2] * v[907] + pM[3][2] * v[908];
		v[450] = pM[3][3] * v[392] + pM[2][3] * v[393] + pM[4][0] * v[841] - pM[0][4] * v[885] + pM[0][3] * v[890]
			+ pM[0][5] * v[892] + pM[4][3] * v[906] - pM[2][4] * v[907] - pM[3][4] * v[908] + pM[2][5] * v[909] + pM[3][5] * v[912]
			+ pM[3][0] * v[924];
		v[414] = pM[1][0] * v[850];
		v[415] = -(pM[0][0] * v[911]);
		v[416] = -(pM[0][4] * v[717]);
		v[417] = pM[1][4] * v[716];
		v[915] = v[414] + v[415] + v[416] + v[417];
		v[518] = pM[1][2] * v[393] - pM[4][4] * v[776] + pM[0][0] * v[817] - pM[4][0] * v[851] + pM[3][0] * v[868]
			- pM[0][2] * v[888] + pM[1][0] * v[897] + pM[4][2] * v[915];
		v[418] = pM[0][3] * v[880];
		v[419] = -(pM[1][3] * v[901]);
		v[913] = v[418] + v[419];
		v[599] = pM[4][0] * v[778] + pM[0][0] * v[807] - pM[1][0] * v[832] + pM[2][0] * (-v[301] + v[845]) + pM[0][2] * v[894]
			- pM[1][2] * v[912] + pM[2][2] * v[913];
		v[593] = -(pM[4][0] * v[779]) - pM[0][0] * v[808] + pM[1][0] * v[831] + pM[3][0] * v[845] - pM[0][2] * v[893]
			+ pM[1][2] * v[909] + pM[3][2] * v[913];
		v[420] = pM[1][0] * v[902];
		v[421] = pM[0][0] * v[914];
		v[916] = v[420] - v[421];
		v[560] = pM[1][0] * v[246] + pM[0][0] * v[292] - pM[4][3] * v[776] + pM[3][0] * v[854] + pM[4][0] * v[855]
			- pM[0][2] * v[887] + pM[1][2] * v[907] + pM[3][2] * v[916] + pM[3][5] * v[922];
		v[469] = pM[1][3] * v[393] + pM[0][4] * v[887] - pM[0][3] * v[888] - pM[0][5] * v[893] - pM[1][4] * v[907]
			+ pM[1][5] * v[909] + pM[3][5] * v[913] + pM[4][3] * v[915] - pM[3][4] * v[916] + pM[4][0] * v[917] + pM[3][0] * v[923];
		v[528] = -(pM[1][2] * v[392]) + pM[4][4] * v[775] + pM[4][0] * v[864] + pM[2][0] * v[868] + pM[0][2] * v[889]
			+ pM[4][2] * v[918] + pM[0][0] * v[919] - pM[1][0] * v[920];
		v[426] = pM[2][3] * v[921];
		v[565] = -(pM[1][0] * v[248]) - pM[4][2] * v[426] + pM[4][3] * v[775] + pM[2][0] * v[854] + pM[4][0] * v[866]
			+ pM[0][2] * v[886] - pM[1][2] * v[908] + pM[2][2] * v[916] + pM[2][5] * v[922];
		v[484] = pM[0][0] * v[340] - pM[1][3] * v[392] + pM[4][4] * v[426] + pM[4][0] * v[862] - pM[0][4] * v[886]
			+ pM[0][3] * v[889] + pM[0][5] * v[894] - pM[1][5] * v[900] + pM[1][4] * v[908] + pM[2][5] * v[913] - pM[2][4] * v[916]
			+ pM[4][3] * v[918] + pM[2][0] * v[923] - pM[1][0] * v[924];
		v[435] = pM[4][0] * v[663];
		v[436] = pM[4][0] * v[666];
		v[438] = pM[1][1] * v[926];
		v[439] = pM[2][1] * v[872];
		v[441] = pM[4][1] * v[651];
		v[932] = -v[435] + v[441] + v[925] - v[928];
		v[442] = pM[4][1] * v[654];
		v[444] = -(pM[3][4] * v[733]);
		v[930] = v[436] - v[438] - v[442] + v[444];
		v[445] = pM[2][0] * v[813];
		v[446] = -(pM[2][4] * v[733]);
		v[931] = -v[439] + v[445] - v[446] + v[927];
		v[582] = -(pM[4][3] * v[789]) + pM[1][1] * v[892] + pM[2][1] * v[893] + pM[3][1] * v[894] + pM[2][3] * v[930]
			+ pM[3][3] * v[931] + pM[1][3] * v[932];
		v[499] = -(pM[3][0] * v[447]) + pM[4][5] * v[929] - pM[2][5] * v[930] - pM[3][5] * v[931] - pM[1][5] * v[932]
			+ pM[1][0] * v[942] + pM[4][4] * (-v[161] + v[962]);
		v[451] = pM[4][0] * v[673];
		v[452] = pM[4][0] * v[667];
		v[547] = pM[2][3] * v[451] + pM[4][3] * v[784] - pM[4][0] * v[842] + pM[0][1] * v[885] + pM[2][1] * v[907]
			+ pM[3][1] * v[908] + pM[3][3] * v[933] + pM[0][3] * v[934] + pM[0][0] * v[935] + pM[0][5] * v[936];
		v[454] = pM[0][1] * v[926];
		v[455] = pM[0][4] * v[937];
		v[456] = pM[0][1] * v[938];
		v[460] = pM[3][4] * v[732];
		v[939] = v[452] - v[454] + v[460];
		v[462] = pM[2][4] * v[732];
		v[940] = -v[455] + v[456] - v[462];
		v[588] = pM[4][3] * (-v[176] + v[791]) + pM[0][1] * v[892] + pM[2][1] * v[909] + pM[3][1] * v[912] + pM[0][3] * v[932]
			+ pM[0][4] * v[936] + pM[2][3] * v[939] + pM[3][3] * v[940];
		v[509] = pM[2][4] * v[451] + pM[4][5] * (-v[170] - v[176] - v[189] - v[650]) - pM[2][0] * v[852] + pM[0][5] * (-v[441]
			+ v[928]) + pM[3][4] * v[933] - pM[2][5] * v[939] - pM[3][5] * v[940] + pM[4][4] * v[941] + pM[0][0] * v[942]
			+ pM[3][0] * v[943];
		v[470] = -(pM[4][0] * v[672]);
		v[472] = pM[4][0] * v[759];
		v[473] = pM[0][4] * v[945];
		v[474] = pM[0][1] * v[872];
		v[475] = -(pM[4][1] * v[661]);
		v[946] = v[451] + v[475];
		v[476] = pM[4][1] * v[655];
		v[970] = -v[476] + v[939];
		v[477] = -(pM[0][5] * v[733]);
		v[478] = -(pM[1][5] * v[732]);
		v[950] = v[472] + v[477] + v[478] + v[944];
		v[561] = pM[0][3] * v[470] + pM[3][0] * v[492] - pM[4][3] * v[785] - pM[4][0] * v[857] - pM[0][1] * v[887]
			+ pM[1][1] * v[907] + pM[3][1] * v[916] + pM[1][3] * v[946] + pM[3][3] * v[950] + pM[3][5] * v[953];
		v[479] = -(pM[0][4] * v[733]);
		v[480] = pM[1][4] * v[732];
		v[947] = -v[473] + v[474] + v[479] - v[480];
		v[594] = -(pM[4][3] * v[792]) - pM[0][1] * v[893] + pM[1][1] * v[909] + pM[3][1] * v[913] - pM[0][3] * v[930]
			+ pM[3][3] * v[947] + pM[1][3] * v[970];
		v[519] = -(pM[1][0] * v[325]) + pM[0][4] * v[470] + pM[1][5] * v[476] + pM[1][4] * v[946] - pM[3][5] * v[947]
			+ pM[4][5] * v[948] - pM[3][0] * v[949] + pM[3][4] * v[950] + pM[4][4] * v[978];
		v[485] = -(pM[4][0] * v[819]);
		v[486] = pM[2][5] * v[945];
		v[954] = v[485] + v[486];
		v[488] = -(pM[0][5] * v[952]);
		v[955] = v[488] - v[933] + v[951];
		v[566] = pM[2][0] * v[492] + pM[4][3] * v[782] + pM[0][1] * v[886] - pM[1][1] * v[908] + pM[2][1] * v[916]
			+ pM[2][3] * v[950] + pM[2][5] * v[953] + pM[0][3] * v[954] + pM[1][3] * v[955];
		v[489] = pM[0][4] * v[952];
		v[983] = v[489] + v[940];
		v[600] = pM[1][3] * (v[455] - v[456] - v[489]) + pM[4][3] * (-v[609] + v[790]) + pM[0][1] * v[894] - pM[1][1] * v[912]
			+ pM[2][1] * v[913] + pM[0][3] * (v[446] + v[931]) + pM[2][3] * v[947] + pM[2][4] * v[953];
		v[529] = pM[4][4] * v[158] - pM[0][0] * v[447] + pM[2][5] * (-v[479] + v[480]) + pM[1][5] * v[489] + pM[4][5] * (-v[616]
			- v[617] - v[636] - v[637]) + pM[1][0] * (-v[463] - v[756]) - pM[2][0] * v[949] + pM[2][4] * v[950] + pM[0][4] * v[954]
			+ pM[1][4] * v[955];
		v[504] = -(pM[3][2] * v[611]);
		v[505] = -(pM[3][2] * v[619]);
		v[960] = v[504] + v[505] + v[956] + v[957] + v[958] + v[959];
		v[584] = -(pM[4][0] * v[800]) - pM[1][2] * v[936] + pM[4][3] * v[960] + pM[1][3] * v[967] + pM[1][0] * v[968]
			+ pM[1][1] * v[969];
		v[583] = -(pM[4][2] * v[789]) + pM[2][2] * v[930] + pM[3][2] * v[931] + pM[1][2] * v[932] + pM[4][4] * v[960];
		v[543] = -(pM[1][1] * v[398]) - pM[2][2] * v[470] + pM[4][0] * (-v[68] - v[82] - v[84]) + pM[1][2] * v[934]
			+ pM[3][2] * v[954] + pM[4][5] * v[960] + pM[4][1] * v[961] + pM[4][2] * (-v[147] + v[962]);
		v[510] = pM[0][2] * v[963];
		v[511] = pM[0][1] * v[964];
		v[512] = pM[0][2] * v[965];
		v[513] = pM[0][0] * v[966];
		v[514] = pM[3][2] * v[610];
		v[515] = pM[3][2] * v[618];
		v[971] = v[510] + v[511] + v[512] + v[513] + v[514] + v[515];
		v[590] = pM[4][0] * v[829] - pM[0][2] * v[936] + pM[0][3] * v[967] + pM[0][0] * v[968] + pM[0][1] * v[969]
			+ pM[4][3] * v[971];
		v[589] = pM[4][2] * v[791] + pM[0][2] * v[932] + pM[2][2] * v[970] + pM[4][4] * v[971] + pM[3][2] * v[983];
		v[548] = -(pM[0][1] * v[398]) + pM[4][1] * (-v[113] - v[686]) + pM[4][2] * (-v[162] + v[784]) + pM[0][2] * v[934]
			+ pM[2][2] * v[946] - pM[3][2] * v[955] + pM[4][5] * v[971] + pM[4][0] * (v[69] + v[972]);
		v[524] = -(pM[3][2] * v[643]);
		v[525] = -(pM[3][2] * v[642]);
		v[977] = v[524] + v[525] + v[973] + v[974] + v[975] + v[976];
		v[596] = -(pM[4][0] * v[846]) + pM[3][1] * v[922] - pM[3][2] * v[953] + pM[4][3] * v[977] + pM[3][3] * v[981]
			+ pM[3][0] * v[982];
		v[595] = -(pM[4][2] * v[792]) - pM[0][2] * v[930] + pM[3][2] * v[947] + pM[1][2] * v[970] + pM[4][4] * v[977];
		v[562] = pM[0][2] * v[470] + pM[3][0] * v[750] + pM[1][2] * v[946] + pM[3][2] * v[950] + pM[4][5] * v[977] + pM[4][2] * (
			-v[659] - v[671] + v[978]) + pM[4][1] * v[979] + pM[4][0] * v[980];
		v[530] = -(pM[0][2] * v[611]);
		v[531] = -(pM[1][2] * v[610]);
		v[532] = -(pM[0][2] * v[619]);
		v[533] = -(pM[1][2] * v[618]);
		v[534] = -(pM[2][2] * v[643]);
		v[535] = -(pM[2][2] * v[642]);
		v[984] = v[530] + v[531] + v[532] + v[533] + v[534] + v[535];
		v[602] = pM[4][0] * v[844] + pM[2][1] * v[922] - pM[2][2] * v[953] + pM[2][3] * v[981] + pM[2][0] * v[982]
			+ pM[4][3] * v[984];
		v[601] = pM[4][2] * v[790] + pM[0][2] * v[931] + pM[2][2] * v[947] - pM[1][2] * v[983] + pM[4][4] * v[984];
		v[567] = pM[4][1] * (-v[612] - v[644] - v[645]) + pM[2][0] * v[750] + pM[4][2] * (-v[156] + v[782]) + pM[2][2] * v[950]
			+ pM[0][2] * v[954] + pM[1][2] * v[955] + pM[4][5] * v[984];
		pMinv[0][0] = v[196] * (pM[5][1] * v[365] - pM[5][2] * v[433] + pM[5][3] * v[497] - pM[5][4] * v[540] + pM[5][5] * v[580]
			);
		pMinv[0][1] = v[196] * (-(pM[5][1] * v[391]) + pM[5][2] * v[449] - pM[5][3] * v[507] + pM[5][4] * v[545]
			- pM[5][5] * v[586]);
		pMinv[0][2] = v[196] * (pM[5][1] * v[412] - pM[5][2] * v[468] + pM[5][3] * v[517] - pM[5][4] * v[559] + pM[5][5] * v[592]
			);
		pMinv[0][3] = v[196] * (-(pM[5][1] * v[423]) + pM[5][2] * v[483] - pM[5][3] * v[527] + pM[5][4] * v[564]
			- pM[5][5] * v[598]);
		pMinv[0][4] = v[196] * (pM[5][1] * v[37] - pM[5][2] * v[62] + pM[5][3] * v[63] - pM[5][4] * v[88] + pM[5][5] * v[89]);
		pMinv[0][5] = -(v[196] * v[363]);
		pMinv[1][0] = v[196] * (-(pM[5][0] * v[365]) + pM[5][2] * v[434] - pM[5][3] * v[498] + pM[5][4] * v[541]
			- pM[5][5] * v[581]);
		pMinv[1][1] = v[196] * (pM[5][0] * v[391] - pM[5][2] * v[450] + pM[5][3] * v[508] - pM[5][4] * v[546] + pM[5][5] * v[587]
			);
		pMinv[1][2] = v[196] * (-(pM[5][0] * v[412]) + pM[5][2] * v[469] - pM[5][3] * v[518] + pM[5][4] * v[560]
			- pM[5][5] * v[593]);
		pMinv[1][3] = v[196] * (pM[5][0] * v[423] - pM[5][2] * v[484] + pM[5][3] * v[528] - pM[5][4] * v[565] + pM[5][5] * v[599]
			);
		pMinv[1][4] = v[196] * (pM[5][4] * v[116] - pM[5][5] * v[141] - pM[5][0] * v[37] + pM[5][2] * v[90] - pM[5][3] * v[91]);
		pMinv[1][5] = v[196] * v[431];
		pMinv[2][0] = v[196] * (pM[5][0] * v[433] - pM[5][1] * v[434] + pM[5][3] * v[499] - pM[5][4] * v[542] + pM[5][5] * v[582]
			);
		pMinv[2][1] = v[196] * (-(pM[5][0] * v[449]) + pM[5][1] * v[450] - pM[5][3] * v[509] + pM[5][4] * v[547]
			- pM[5][5] * v[588]);
		pMinv[2][2] = v[196] * (pM[5][0] * v[468] - pM[5][1] * v[469] + pM[5][3] * v[519] - pM[5][4] * v[561] + pM[5][5] * v[594]
			);
		pMinv[2][3] = v[196] * (-(pM[5][0] * v[483]) + pM[5][1] * v[484] - pM[5][3] * v[529] + pM[5][4] * v[566]
			- pM[5][5] * v[600]);
		pMinv[2][4] = v[196] * (pM[5][3] * v[142] - pM[5][4] * v[167] + pM[5][5] * v[168] + pM[5][0] * v[62] - pM[5][1] * v[90]);
		pMinv[2][5] = -(v[196] * v[495]);
		pMinv[3][0] = v[196] * (-(pM[5][0] * v[497]) + pM[5][1] * v[498] - pM[5][2] * v[499] + pM[5][4] * v[543]
			- pM[5][5] * v[583]);
		pMinv[3][1] = v[196] * (pM[5][0] * v[507] - pM[5][1] * v[508] + pM[5][2] * v[509] - pM[5][4] * v[548] + pM[5][5] * v[589]
			);
		pMinv[3][2] = v[196] * (-(pM[5][0] * v[517]) + pM[5][1] * v[518] - pM[5][2] * v[519] + pM[5][4] * v[562]
			- pM[5][5] * v[595]);
		pMinv[3][3] = v[196] * (pM[5][0] * v[527] - pM[5][1] * v[528] + pM[5][2] * v[529] - pM[5][4] * v[567] + pM[5][5] * v[601]
			);
		pMinv[3][4] = v[196] * (-(pM[5][2] * v[142]) + pM[5][4] * v[169] - pM[5][5] * v[194] - pM[5][0] * v[63]
			+ pM[5][1] * v[91]);
		pMinv[3][5] = v[196] * v[538];
		pMinv[4][0] = v[196] * (pM[5][0] * v[540] - pM[5][1] * v[541] + pM[5][2] * v[542] - pM[5][3] * v[543] + pM[5][5] * v[584]
			);
		pMinv[4][1] = v[196] * (-(pM[5][0] * v[545]) + pM[5][1] * v[546] - pM[5][2] * v[547] + pM[5][3] * v[548]
			- pM[5][5] * v[590]);
		pMinv[4][2] = v[196] * (pM[5][0] * v[559] - pM[5][1] * v[560] + pM[5][2] * v[561] - pM[5][3] * v[562] + pM[5][5] * v[596]
			);
		pMinv[4][3] = v[196] * (-(pM[5][0] * v[564]) + pM[5][1] * v[565] - pM[5][2] * v[566] + pM[5][3] * v[567]
			- pM[5][5] * v[602]);
		pMinv[4][4] = v[196] * (-(pM[5][1] * v[116]) + pM[5][2] * v[167] - pM[5][3] * v[169] + pM[5][5] * v[195]
			+ pM[5][0] * v[88]);
		pMinv[4][5] = -(v[196] * v[578]);
		pMinv[5][0] = v[196] * (-(pM[5][0] * v[580]) + pM[5][1] * v[581] - pM[5][2] * v[582] + pM[5][3] * v[583]
			- pM[5][4] * v[584]);
		pMinv[5][1] = v[196] * (pM[5][0] * v[586] - pM[5][1] * v[587] + pM[5][2] * v[588] - pM[5][3] * v[589] + pM[5][4] * v[590]
			);
		pMinv[5][2] = v[196] * (-(pM[5][0] * v[592]) + pM[5][1] * v[593] - pM[5][2] * v[594] + pM[5][3] * v[595]
			- pM[5][4] * v[596]);
		pMinv[5][3] = v[196] * (pM[5][0] * v[598] - pM[5][1] * v[599] + pM[5][2] * v[600] - pM[5][3] * v[601] + pM[5][4] * v[602]
			);
		pMinv[5][4] = v[196] * (pM[5][1] * v[141] - pM[5][2] * v[168] + pM[5][3] * v[194] - pM[5][4] * v[195] - pM[5][0] * v[89]
			);
		pMinv[5][5] = v[196] * v[605];
#pragma endregion

		Minv.PtrToMatrix(pMinv, 6);

		for (int i = 0; i < 6; i++)
			delete[] pM[i];
		delete[] pM;

		for (int i = 0; i < 6; i++)
			delete[] pMinv[i];
		delete[] pMinv;

		return Minv;
	}
}

//Operador produto tensorial entre dois vetores
Matrix dyadic(Matrix &matrix1, Matrix &matrix2)
{
	//Verificação da possibilidade do produto
	if (matrix1.m_columns != 1 || matrix2.m_columns != 1 || matrix1.m_lines != matrix2.m_lines)
	{
		printf("Nao e possivel calcular o produto tensorial. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		int order = matrix1.m_lines;
		Matrix return_m(order, order);
		for (int i = 0; i < order; i++)
		for (int j = 0; j < order; j++)
				return_m(i,j)=matrix1(i,0)*matrix2(j,0);
		return return_m;
	}
}
//Operador skew de um vetor
Matrix skew(Matrix &matrix1)
{
	//Verificação da possibilidade do produto
	if (matrix1.m_columns != 1 || matrix1.m_lines != 3)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		Matrix return_m(3,3);
		return_m(0,1)=-matrix1(2,0);
		return_m(0,2)=+matrix1(1,0);
		return_m(1,2)=-matrix1(0,0);

		return_m(1,0)=+matrix1(2,0);
		return_m(2,0)=-matrix1(1,0);
		return_m(2,1)=+matrix1(0,0);
		return return_m;
	}
}
//Operador axial de um vetor
Matrix axial(Matrix &matrix1)
{
	//Verificação da possibilidade do produto
	if (matrix1.m_columns != 3 || matrix1.m_lines != 3)
	{
		printf("Nao e possivel calcular o produto escalar. Dimensoes incompativeis!");
		return 0;
	}
	else
	{
		Matrix return_m(3);
		return_m(0,0)=-matrix1(1,2);
		return_m(1,0)=+matrix1(0,2);
		return_m(2,0)=-matrix1(0,1);
		return return_m;
	}
}
//Resolve o sistema linear da forma Ax=b (utilizando fatoração LU)
Matrix fullsystem(Matrix &A, Matrix &b, int *flag_error)
{
	int n = A.getLines();
	int nrhs = 1;
	int lda = A.getLines();
	int ldb = b.getLines();
	int info = 1;
	int *ipiv;
	char param = 'N';
	ipiv = new int[b.getLines()];
	//Faz fatoracao LU
	dgetrf(&n,&n,A.getMatrix(),&lda,ipiv,&info);
	if (info != 0)
		*flag_error = 1;
	else
		*flag_error = 0;
	//Resolve o sistema linear
	dgetrs(&param,&n,&nrhs,A.getMatrix(),&lda,ipiv,b.m_matrix,&ldb,&info);
	if (info != 0)
		*flag_error = 1;
	else
	{
		if (*flag_error == 0)
			*flag_error = 0;
	}
		
	delete []ipiv;
	return b;
}
//Retorna o endereço de uma matriz
double* Matrix::getMatrix()
{
	return m_matrix;
}

//Retorna a norma de um vetor
double norm(Matrix &matrix1)
{
	if (matrix1.getColumns() != 1)
		printf("Dimensao nao consistente para calculo da norma");
	else
	{
		if (matrix1.getLines() != 3)
		{
			//Marina
			if (matrix1.getLines() == 1)
			{
				double return_value = sqrt(matrix1(0, 0)*matrix1(0, 0));
				return return_value;
			}
			if (matrix1.getLines() == 2)
			{
				double return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
					matrix1(1, 0)*matrix1(1, 0));
				return return_value;
			}
			if (matrix1.getLines() == 4)
			{
				double return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
					matrix1(1, 0)*matrix1(1, 0) +
					matrix1(2, 0)*matrix1(2, 0) + 
					matrix1(3, 0)*matrix1(3, 0) );
				return return_value;
			}
			if (matrix1.getLines() == 6)
			{
				double return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
					matrix1(1, 0)*matrix1(1, 0) +
					matrix1(2, 0)*matrix1(2, 0) +
					matrix1(3, 0)*matrix1(3, 0) + 
					matrix1(4, 0)*matrix1(4, 0) +
					matrix1(5, 0)*matrix1(5, 0));
				return return_value;
			}
			//printf("Norma infinito\n");
			double max = 0;
			for (int i=0; i< matrix1.getLines(); i++)
			{
				//Detecção de NaN
				if (matrix1(i, 0) == matrix1(i, 0))
				{
					//Detecção de infinito
					if (matrix1(i, 0) >= 1e300 || matrix1(i, 0) <= -1e300)
						return 1e100;
					else
					{
						if (max < abs(matrix1(i, 0)))
						{
							max = abs(matrix1(i, 0));
							//printf("GL %d\n",i + 1);
						}
							
					}
					
				}
				else
					return 1e100;//valor muito alto, pois detectou NaN
				
			}
			return max;
		}
		else
		{
			double return_value = sqrt( matrix1(0,0)*matrix1(0,0) + 
										matrix1(1,0)*matrix1(1,0) +
										matrix1(2,0)*matrix1(2,0) );
			return return_value;
		}
	}
	return 0;
}

//Retorna a norma de um vetor considerando somente os 4 primeiros graus de liberdade
double norm4(Matrix &matrix1)
{
	if (matrix1.getLines() < 4)
		printf("Error. Function norm4\n");
	else
	{
		double return_value = sqrt(matrix1(0, 0)*matrix1(0, 0) +
			matrix1(1, 0)*matrix1(1, 0) +
			matrix1(2, 0)*matrix1(2, 0) +
			matrix1(3, 0)*matrix1(3, 0));
		return return_value;
	}
	return 0;
}

//Retorna a transposta de uma matriz
Matrix transp(Matrix &matrix1)
{
	Matrix answer(matrix1.getColumns(),matrix1.getLines());
	for (int j=0;j<matrix1.getColumns();j++)
	{
		for (int i=0;i<matrix1.getLines();i++)
		{	
			answer(j,i)=matrix1(i,j);
		}
	}
	return answer;
}
//Zera a matriz
void zeros(Matrix* matrix1)
{
	for (int j=0;j<matrix1->getColumns();j++)
	{
		for (int i=0;i<matrix1->getLines();i++)
		{	
			(*matrix1)(i,j) = 0.0;
		}
	}
}
//Função para o calculo do operador V (para montagem da matriz de rigidez geometrica)
Matrix V(Matrix x, Matrix t,double alpha_escalar)
{
	//Funçoes da pag. 50 e afins da tese M. Lourdes
	double h = 4.0/(4.0+alpha_escalar*alpha_escalar);
	double h2 = 0.5*h;
	double h3=0;
	double h4=-0.25*h*h;
	double h5=0;
	double h8=-0.5*h*h;

	Matrix return_m = dyadic(h8*t-h4*(skew(x)*t)+h5*(skew(x)*skew(x))*t,x) + 
		 h2*skew(t)-h3*(2*skew(x)*skew(t) - skew(t)*skew(x));
	return return_m;
}
//Função para o calculo do operador d_V (para montagem da matriz de rigidez geometrica)
Matrix d_V(Matrix x,Matrix d_x, Matrix t,double alpha_escalar)
{
	//Funçoes da pag. 50 e afins da tese M. Lourdes
	double h = 4.0/(4.0+alpha_escalar*alpha_escalar);
	double h3=0;
	double h4=-0.25*h*h;
	double h5=0;
	double h6=0.25*h*h*h;
	double h7=0;
	double h8=-0.5*h*h;
	double h9=0.5*h*h*h;

	Matrix return_m = dot(x,d_x)*(dyadic(h9*t-h6*(skew(x)*t)+h7*(skew(x)*skew(x))*t,x)) +
		dyadic(h8*t-h4*skew(x)*t+h5*(skew(x)*skew(x))*t,d_x) +
		dyadic(h5*(skew(x)*skew(d_x) + skew(d_x)*skew(x))*t - h4*skew(d_x)*t,x) +
		h4*dot(x,d_x)*skew(t) - h5*dot(x,d_x)*(2*skew(x)*skew(t)-skew(t)*skew(x)) -
		h3*(2*skew(d_x)*skew(t) - skew(t)*skew(d_x));
	return return_m;
}

//Retorna um vetor com esses componentes a,b,c
Matrix List(double a, double b, double c)
{
	Matrix return_m(3);
	return_m(0, 0) = a;
	return_m(1, 0) = b;
	return_m(2, 0) = c;
	return return_m;
}
double Power(double a, double b)
{
	return pow(a, b);
}
double Power(Matrix a, double b)
{
	if (a.m_lines != 1 || a.m_columns != 1)
	{
		if (a.m_lines == 3 && a.m_columns == 1)
		{
			if (b == 2)
				return dot(a,a);
			else
				printf("Error in Power function. Supposed to receive a 1x1 or 3x1 matrix!\n");
		}	
		else
			printf("Error in Power function. Supposed to receive a 1x1 or 3x1 matrix!\n");
	}

	return pow(a(0,0), b);
}
double Sin(double a)
{
	return sin(a);
}
double Cos(double a)
{
	return cos(a);
}
//Operador Soma
double operator + (double a, Matrix &matrix2)
{
	//Verifica se as dimensões das matrizes são compativeis
	if ((matrix2.getLines() != 1) || (matrix2.getColumns() != 1))
	{
		//Mensagem de erro - Dimensões Incompativeis
		printf("Matriz deve ser unitaria! Função operator+(double,matrix)\n");
		//Retorna vazio
		return NULL;
	}
	//Caso as dimensões sejam compativeis
	else
	{
		return matrix2(0,0)+a;
	}
}

//Produto de matrizes
Matrix Dot(Matrix &matrix1, Matrix &matrix2)
{
	return matrix1*matrix2;
}

//Retorna um ponteiro de duas dimensões que aponta para a matriz  
void Matrix::MatrixToPtr(double** ptr, int order)
{
	for (int i = 0; i < order; i++)
	{
		for (int j = 0; j < order; j++)
			ptr[i][j] = m_matrix[i + j*order];
	}
	
}

//Salva na matrix o conteudo do double**
void Matrix::PtrToMatrix(double** ptr, int order)
{
	for (int i = 0; i < order; i++)
	{
		for (int j = 0; j < order; j++)
			m_matrix[i + j*order] = ptr[i][j];
	}
	
}

//Salva na matrix o conteudo do double**
void Matrix::PtrToMatrix(double** ptr, int lines, int columns)
{
	for (int i = 0; i < lines; i++)
	{
		for (int j = 0; j < columns; j++)
			m_matrix[i + j*lines] = ptr[i][j];
	}
}



//Calcula os autovalores e autovetores de uma matriz simetrica
//A: matriz de entrada
//P: matriz de saida - possui autovetores organizados em colunas
//D: matrix de saida - possui autovalores na diagonal principal
//abstol: tolerancia
int fulleigen1(Matrix &A, Matrix &P, Matrix &D, double abstol)
{
	char jobz = 'V';
	char uplo = 'U';
	char range = 'A';
	double vl, vu;
	int il, iu;
	int m;
	int n = A.getLines();
	Matrix d(n);//onde autovalores são salvos
	int lwork = n * 26;
	Matrix work(lwork);
	int liwork = 10 * n;
	int *iwork;
	iwork = new int[liwork];
	int info;
	Matrix z(n, n);
	int lisuppz = 2 * n;
	int *isuppz;
	isuppz = new int[lisuppz];
	dsyevr(&jobz, &range, &uplo, &n, A.getMatrix(), &n, &vl, &vu, &il, &iu, &abstol, &m, d.getMatrix(), z.getMatrix(), &n, isuppz,work.getMatrix(),&lwork,iwork,&liwork,&info);
	
	if (info != 0)
	{
		printf("Error in function 'fulleigen'!\n");
		return info;
	}
	//Valores para retorno:
	for (int i = 0; i < n; i++)
	{
		//Autovalores:
		D(i, i) = d(i, 0);
		//Autovetores ja estão salvos na matriz P
	}
	P = z;//Copia valores da saida (autovetores) para a matriz P - organizados em colunas

	delete []isuppz;
	delete []iwork;
	return info;
}
//Calcula os autovalores e autovetores de uma matriz simetrica
//A: matriz de entrada
//P: matriz de saida - possui autovetores organizados em colunas
//d: matrix de saida - possui autovalores na diagonal principal
int fulleigen2(Matrix &A, Matrix &P, Matrix &D)
{
	char jobz = 'V';
	char uplo = 'U';
	int n = A.getLines();
	Matrix d(n);//onde autovalores são salvos
	int lwork = n * 3;
	Matrix work(lwork);
	int info;
	dsyev(&jobz, &uplo, &n, A.getMatrix(), &n, d.getMatrix(), work.getMatrix(), &lwork, &info);

	if (info != 0)
	{
		printf("Error in function 'eigen'!\n");
		return info;
	}
	//Valores para retorno:
	for (int i = 0; i < n; i++)
	{
		//Autovalores:
		D(i, i) = d(i, 0);
		//Autovetores ja estão salvos na matriz P
	}
	P = A;//Copia valores da saida (autovetores) para a matriz P - organizados em colunas
	return info;
}
//Calcula o menor autovalor de uma matriz simetrica
double mineigen(Matrix &A, Matrix &P, Matrix &D, double abstol)
{
	char jobz = 'N';
	char uplo = 'U';
	char range = 'I';
	double vl, vu;
	int il = 1;
	int iu = 1;
	int m;
	int n = A.getLines();
	Matrix d(n);//onde autovalores são salvos
	int lwork = n * 26;
	Matrix work(lwork);
	int liwork = 10 * n;
	int *iwork;
	iwork = new int[liwork];
	int info;
	Matrix z(1);
	int lisuppz = 2 * n;
	int *isuppz;
	isuppz = new int[lisuppz];
	dsyevr(&jobz, &range, &uplo, &n, A.getMatrix(), &n, &vl, &vu, &il, &iu, &abstol, &m, d.getMatrix(), z.getMatrix(), &n, isuppz, work.getMatrix(), &lwork, iwork, &liwork, &info);

	if (info != 0)
	{
		printf("Error in function 'fulleigen'!\n");
	}
	//Valores para retorno:
	for (int i = 0; i < n; i++)
	{
		//Autovalores:
		D(i, i) = d(i, 0);
		//Autovetores ja estão salvos na matriz P
	}
	P = z;//Copia valores da saida (autovetores) para a matriz P - organizados em colunas

	delete[]isuppz;
	delete[]iwork;
	return d(0,0);
}

//Calcula redução a primeira volta [-pi,pi]
double ArcReduction(double arc)
{
	double cosarc = cos(arc);
	double sinarc = sin(arc);

	//Evita erros numericos
	if (cosarc > 1.0)
		cosarc = 1.0;
	if (cosarc < -1.0)
		cosarc = -1.0;
	if (sinarc > 1.0)
		sinarc = 1.0;
	if (sinarc < -1.0)
		sinarc = -1.0;

	double ret_arc;

	if (sinarc >= 0 && cosarc >= 0)
		ret_arc = asin(sinarc);
	if (sinarc >= 0 && cosarc < 0)
		ret_arc = -asin(sinarc) + PI;
	if (sinarc < 0 && cosarc >= 0)
		ret_arc = asin(sinarc);
	if (sinarc < 0 && cosarc < 0)
		ret_arc = -asin(sinarc) - PI;

	return ret_arc;
}

//Calcula redução a primeira volta [0,2*pi]
double ArcReduction2p(double arc)
{
	double cosarc = cos(arc);
	double sinarc = sin(arc);

	//Evita erros numericos
	if (cosarc > 1.0)
		cosarc = 1.0;
	if (cosarc < -1.0)
		cosarc = -1.0;
	if (sinarc > 1.0)
		sinarc = 1.0;
	if (sinarc < -1.0)
		sinarc = -1.0;

	double ret_arc;

	if (sinarc >= 0 && cosarc >= 0)
		ret_arc = asin(sinarc);
	if (sinarc >= 0 && cosarc < 0)
		ret_arc = -asin(sinarc) + PI;
	if (sinarc < 0 && cosarc >= 0)
		ret_arc = asin(sinarc) + 2 * PI;
	if (sinarc < 0 && cosarc < 0)
		ret_arc = -asin(sinarc) - PI + 2 * PI;

	return ret_arc;
}

//Marina
//Resolve polinomio de 3o grau
 //---------------------------------------------------------------------------
 // solve cubic equation x^3 + a*x^2 + b*x + c
 // x - array of size 3
 // In case 3 real roots: => x[0], x[1], x[2], return 3
 //         2 real roots: x[0], x[1],          return 2
 //         1 real root : x[0], x[1] ± i*x[2], return 1
unsigned int solveP3(double *x, double a, double b, double c) {
	double a2 = a * a;
	double q = (a2 - 3 * b) / 9;
	double r = (a*(2 * a2 - 9 * b) + 27 * c) / 54;
	double r2 = r * r;
	double q3 = q * q*q;
	double A, B;
	if (r2 < q3)
	{
		double t = r / sqrt(q3);
		if (t < -1) t = -1;
		if (t > 1) t = 1;
		t = acos(t);
		a /= 3; q = -2 * sqrt(q);
		x[0] = q * cos(t / 3) - a;
		x[1] = q * cos((t + M_2PI) / 3) - a;
		x[2] = q * cos((t - M_2PI) / 3) - a;
		return 3;
	}
	else
	{
		A = -pow(fabs(r) + sqrt(r2 - q3), 1. / 3);
		if (r < 0) A = -A;
		B = (0 == A ? 0 : q / A);

		a /= 3;
		x[0] = (A + B) - a;
		x[1] = -0.5*(A + B) - a;
		x[2] = 0.5*sqrt(3.)*(A - B);
		if (fabs(x[2]) < eps) { x[2] = x[1]; return 2; }

		return 1;
	}
}

//Resolve polinômio de 4o grau
//---------------------------------------------------------------------------
// solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d
// Attention - this function returns dynamically allocated array. It has to be released afterwards.
DComplex* solve_quartic(double a, double b, double c, double d)
{
	double a3 = -b;
	double b3 = a * c - 4.*d;
	double c3 = -a * a*d - c * c + 4.*b*d;

	// cubic resolvent
	// y^3 - b*y^2 + (ac-4d)*y - a^2*d-c^2+4*b*d = 0

	double x3[3];
	unsigned int iZeroes = solveP3(x3, a3, b3, c3);

	double q1, q2, p1, p2, D, sqD, y;

	y = x3[0];
	// The essence - choosing Y with maximal absolute value.
	if (iZeroes != 1)
	{
		if (fabs(x3[1]) > fabs(y)) y = x3[1];
		if (fabs(x3[2]) > fabs(y)) y = x3[2];
	}

	// h1+h2 = y && h1*h2 = d  <=>  h^2 -y*h + d = 0    (h === q)

	D = y * y - 4 * d;
	if (fabs(D) < eps) //in other words - D==0
	{
		q1 = q2 = y * 0.5;
		// g1+g2 = a && g1+g2 = b-y   <=>   g^2 - a*g + b-y = 0    (p === g)
		D = a * a - 4 * (b - y);
		if (fabs(D) < eps) //in other words - D==0
			p1 = p2 = a * 0.5;

		else
		{
			sqD = sqrt(D);
			p1 = (a + sqD) * 0.5;
			p2 = (a - sqD) * 0.5;
		}
	}
	else
	{
		sqD = sqrt(D);
		q1 = (y + sqD) * 0.5;
		q2 = (y - sqD) * 0.5;
		// g1+g2 = a && g1*h2 + g2*h1 = c       ( && g === p )  Krammer
		p1 = (a*q1 - c) / (q1 - q2);
		p2 = (c - a * q2) / (q1 - q2);
	}

	DComplex* retval = new DComplex[4];

	// solving quadratic eq. - x^2 + p1*x + q1 = 0
	D = p1 * p1 - 4 * q1;
	if (D < 0.0)
	{
		retval[0].real(-p1 * 0.5);
		retval[0].imag(sqrt(-D) * 0.5);
		retval[1] = std::conj(retval[0]);
	}
	else
	{
		sqD = sqrt(D);
		retval[0].real((-p1 + sqD) * 0.5);
		retval[1].real((-p1 - sqD) * 0.5);
	}

	// solving quadratic eq. - x^2 + p2*x + q2 = 0
	D = p2 * p2 - 4 * q2;
	if (D < 0.0)
	{
		retval[2].real(-p2 * 0.5);
		retval[2].imag(sqrt(-D) * 0.5);
		retval[3] = std::conj(retval[2]);
	}
	else
	{
		sqD = sqrt(D);
		retval[2].real((-p2 + sqD) * 0.5);
		retval[3].real((-p2 - sqD) * 0.5);
	}

	return retval;
}
