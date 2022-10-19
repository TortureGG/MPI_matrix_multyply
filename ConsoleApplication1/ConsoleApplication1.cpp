#include <iostream>
#include "mpi.h"
#include "math.h"
using namespace std;

class Table {
	public:
		double** Matrix;
		int wide, height;

        int Q, S;
        int* Compose_Q;
        int* Compose_S;
        

        void init(int w, int h) {

            this->wide = w;
            this->height = h;

            Matrix = new double* [h];
            for (int i = 0; i < h; i++) {
                Matrix[i] = new double[w];
            }
        }

        void create(int w, int h, int * Compose_Q, int * Compose_S) {

            init(w, h);

            this->Compose_Q = Compose_Q;
            this->Compose_S = Compose_S;


			int count = 0;
			for (int i = 0; i < this->height; i++) {
				for (int j = 0; j < this->wide; j++) {
                    Matrix[i][j] = count;
					count++;
				}
			}
		}

        void fillTableForProc(int* Compose_Q, int* Compose_S, int Q, int S, int id, bool Reversed) {

            this->Q = Q;
            this->S = S;
            this->Compose_Q = Compose_Q;
            this->Compose_S = Compose_S;

            //приплюсовываем строки для нижнего ряда процессоров
            //для этого вычислим целую часть процессора от ширины 
            int num = 0;
            int cel = id / this->S;
            int idEnd = id - cel * this->S;

            if (Reversed) {
                this->Q = S;
                this->S = Q;
                this->Compose_Q = Compose_S;
                this->Compose_S = Compose_Q;
                cel = id - cel * S;
                idEnd = id / S;
            }

            //printf("Celaya chast = %i\n", cel);

            for (int i = 0; i < cel; i++) {
                for (int l = 0; l < S; l++) {
                    num += this->Compose_Q[i] * this->Compose_S[l];
                }
            }
            //printf("num = %i\n", num);

            for (int i = 0; i < this->height; i++) {

                //прибавляем если слева и справа отступы которые забрали другие процессоры
                for (int l = 0; l < idEnd; l++) {
                    num += this->Compose_S[l];
                }

                for (int j = 0; j < this->wide; j++) {
                    Matrix[i][j] = num;
                    num++;
                }
                
                for (int l = idEnd + 1; l < S; l++) {
                    num += this->Compose_S[l];
                }
            }
        }
        
        //необязательно, так как выше сделал Reversed
        void fillScaleTableForProc(int* Compose_Q, int* Compose_S, int Q, int S, int id) {
            //подобное только idEnd и cel поменять местами
            this->Compose_Q = Compose_Q;
            this->Compose_S = Compose_S;
            
            int num = 0;
            int cel = id / S;
            int idEnd = id - cel * S;

            for (int i = 0; i < idEnd; i++) {
                for (int l = 0; l < Q; l++) {
                    num += this->Compose_S[i] * this->Compose_Q[l];
                }
            }

            for (int i = 0; i < this->height; i++) {

                for (int l = 0; l < cel; l++) {
                    num += this->Compose_Q[l];
                }

                for (int j = 0; j < this->wide; j++) {
                    Matrix[i][j] = num;
                    num++;
                }

                for (int l = cel + 1; l < Q; l++) {
                    num += this->Compose_Q[l];
                }
            }

  
        }

        void print() {
            for (int i = 0; i < this->height; i++) {
                for (int j = 0; j < this->wide; j++) {
                    printf("%f\t", Matrix[i][j]);
                }
                printf("\n");
            }
            printf("\n");
        }
		void printForMaintable() {

            int comp[2]  = { 0, 0 };
            int step[2]  = { Compose_Q[comp[0]], Compose_S[comp[1]] };
            int count[2] = { 0, 0 };

            string stuff(19* this->wide, '=');
            cout << "\n" << stuff << endl;

			for (int i = 0; i < this->height; i++) {
                count[0]++;
                count[1] = 0;
                comp[1]  = 0;
                step[1]  = Compose_S[comp[1]];

				for (int j = 0; j < this->wide; j++) {
                    count[1]++;
					printf("%f\t", Matrix[i][j]);
                    if (count[1] == step[1]) {
                        printf(" ||\t");
                        comp[1]++;
                        count[1] = 0;
                        step[1] = Compose_S[comp[1]];
                    }
				}
                if (count[0] == step[0]) {
                    cout << "\n" << stuff << endl;
                    comp[0]++;
                    count[0] = 0;
                    step[0] = Compose_Q[comp[0]];
                }

				printf("\n");
			}
		}

		void clear() {
			for (int i = 0; i < this->height; i++)
				delete[] Matrix[i];

			delete[] Matrix;
		}
};

class Config {

public:
	int id, numproc;
    int wide, height;
    int wideProc, heightProc;
    int lustProc;

	int Q, S;
    int* Compose_Q = new int[Q];
    int* Compose_S = new int[S];

	void set(int numproc, int ID, int wide, int height) {
        this->id = ID;
        this->numproc = numproc;
        this->wide = wide;
        this->height = height;
		//найдем разбиения отталкиваясь от кол-ва процессоров: 6 = 2*3 или 12 = 3*4

        //if (id == 0) printf("Choose Q and S composition:\n");
		int max = 99999;
		for (int i = 1; i <= this->numproc; ++i)
			if (this->numproc % i == 0 && (i + this->numproc / i) < max) {
				Q = i;
				S = this->numproc / i;
				max = Q + S;
                //if (id == 0) printf("N=%d %d * %d = %d \n", this->numproc, Q, S, max);
			}
        //Если у нас композиция вышла больше чем строк и столбцов либо очень много процессоров
        if (Q > this->height) {
            printf("N=%d %d * %d = %d \n", this->numproc, Q, S, max);
            printf("Q > height = %i > %i\n", Q, this->height);
            set(numproc-1, ID, wide, height);
            return;
        }
        if (S > this->wide) {
            printf("N=%d %d * %d = %d \n", this->numproc, Q, S, max);
            printf("S > wide = %i > %i\n", S, this->wide);
            set(numproc-1, ID, wide, height);
            return;
        }
		//после нужно найти шаг разбиения для процессоров
        Compose_Q = decompose(this->height, Q);
        Compose_S = decompose(this->wide, S);

        printf("Choose Q and S composition:\n");
        printf("N=%d %d * %d = %d \n", this->numproc, Q, S, max);

        if (id == 0) {
           /* printf("Choose Q and S composition:\n");
            printf("N=%d %d * %d = %d \n", this->numproc, Q, S, max);*/

            printf("\nComposition steps in main table:\n");
            printf("\nCompose_Q: "); for (int i = 0; i < Q; i++) printf("%d\t", Compose_Q[i]);
            printf("\nCompose_S: "); for (int i = 0; i < S; i++) printf("%d\t", Compose_S[i]);
            printf("\n");
        }
        

        {
            printf("\nTable Height and Wide for each proc:\n");
            int id = 0;
            for (int i = 0; i < Q; i++) {
                for (int j = 0; j < S; j++) {
                    if (this->id == id) {
                        printf("[%i]: %d x %d", id, Compose_Q[i], Compose_S[j]);
                        this->heightProc = Compose_Q[i];
                        this->wideProc = Compose_S[j];
                    }
                    id++;
                }
            }
            printf("\n");
        }

	}

    int* decompose(int Chislo, int amount) {


        int* Compose = new int[Chislo];
        int* Result = new int[amount];

        double val = (double)Chislo / (double)amount;
        double fractpart, intpart;

        fractpart = modf(val, &intpart);
        //printf("%f %f %f", val, intpart, fractpart);

        if (fractpart == 0.0) {
            for (int i = 0; i < amount; i++) {
                Result[i] = val;
            }

            return Result;
        }

        for (int i = 0; i < Chislo; i++) {
            Compose[i] = 1;
            //printf("%d\t", Compose[i]);
        }
        //printf("\n");

        int size = Chislo;
        int num = 1;

        while (size > amount) {

            for (int i = 0; i < Chislo; i++) {
                if (size <= amount) break;
                //if (Compose[i] == 0) continue;
                //ищем наименьшее для сложения
                if (Compose[i] == num) {
                    //ищем последующее похожее
                    for (int j = 0; j < Chislo; j++) {
                        //printf("\nI=%d J=%d S=%d N=%d ",i, j, size, num);
                        /*if (i == j) continue;*/

                        if (Compose[j] == num && i != j) {

                            Compose[i] = Compose[i] + Compose[j];
                            Compose[j] = 0;
                            size--;
                            /*printf("Combine: ");
                            for (int i = 0; i < Chislo; i++) {
                                printf("%d\t", Compose[i]);
                            }
                            printf("\n");*/
                            break;
                        }
                        else if (j == Chislo - 1) {
                            //не нашли
                            num++;
                            j = -1;
                            //printf("\nNot faund J=0\n");
                            continue;
                        }
                    };
                }
                else if (i == Chislo - 1 && size > 1) {
                    num++;
                    i = -1;
                    continue;
                }
            }
        }

        //перекидываем нули в конец
        bool change = true;
        while (change) {
            change = false;
            for (int i = 0; i < Chislo - 1; i++) {

                if (Compose[i] == 0 && Compose[i + 1] != 0) {
                    int buf = Compose[i];
                    Compose[i] = Compose[i + 1];
                    Compose[i + 1] = buf;
                    change = true;
                    //break;
                }
            }
        }

        int count = 0;
        //printf("Sorted: ");
        //for (int i = 0; i < Chislo; i++) {
        //    if (Compose[i] != 0) count++;
        //    printf("%d\t", Compose[i]);
        //}
        //printf(" = %d\n", count);

        //printf("Result: ");
        for (int i = 0; i < amount; i++) {
            Result[i] = Compose[i];
            //printf("%d\t", Result[i]);
        }
        //printf("\n");

        return Result;
    }

};

class Proc{
public:
    int id;
    Table MainTable;
    Table ScaleTable;
    double** AllRowRecv;
    int* neighbourID;
    double* sum;//сумма которая отправляется в последний процеессор, нужна для вывода в последнем процессоре

    void set(int ID, int wide, int height) {
        this->id = ID;

        //таблицы либо получить либо заполнить на основе разбиений
        MainTable.init(wide, height);
        ScaleTable.init(height, wide);

        for (int i = 0; i < height; i++) {
            for (int j = 0; j < wide; j++) {
                MainTable.Matrix[i][j] = 0;
                ScaleTable.Matrix[j][i] = 0;
            }
        }
    }

    void fillTable(int* Compose_Q, int* Compose_S, int Q, int S) {
        MainTable.fillTableForProc(Compose_Q, Compose_S, Q, S, this->id, false);
        //ScaleTable.fillTableForProc(Compose_Q, Compose_S, Q, S, this->id, true);
        ScaleTable.fillScaleTableForProc(Compose_Q, Compose_S, Q, S, this->id);
    }

    void calculateNeighbourID(int Q, int S) {
        int cel = this->id / S;
        int startProc = this->id - cel * S;
        neighbourID = new int[Q];

        for (int i = 0; i < Q; i++) {
            neighbourID[i] = startProc;
            startProc += S;
        }

        //вывод каким процессорам отправляет
        {
            printf("Neighbour Proc: ");
            for (int i = 0; i < Q; i++) {
                if (this->id == neighbourID[i]) printf("[");
                                                printf("%i", neighbourID[i]);
                if (this->id == neighbourID[i]) printf("]");
                printf("\t");
            }
            printf("\n");
        }
    }

    void sendScaleRow(int numproc, int* Compose_Q, int* Compose_S, int Q, int S) {
        //////////////////////////////////////////////////////////////////////////////////////////
        int lenSend = ScaleTable.wide * ScaleTable.height;
        double* RowSend = new double[lenSend];
        int iter = 0;
        for (int j = 0; j < ScaleTable.wide; j++) {
            for (int i = 0; i < ScaleTable.height; i++) {
                RowSend[iter] = ScaleTable.Matrix[i][j];
                iter++;
            }
        }
        printf("Row: "); for (int i = 0; i < lenSend; i++) printf("%f ", RowSend[i]); printf("\n");
        //////////////////////////////////////////////////////////////////////////////////////////
        //рассылка
        for (int i = 0; i < Q; i++) {
            send(RowSend, lenSend, neighbourID[i], 0);
            //if (this->id != neighbourID[i]) {
            //    send(RowSend, lenSend, neighbourID[i]);
            //};
        }
    }
    
    void recvScaleRow(int numproc, int* Compose_Q, int* Compose_S, int Q, int S){
        AllRowRecv = new double* [Q];

        //получаем
        for (int i = 0; i < Q; i++) {

            ///////////////////////////////////////////////////////////////////////
            //найдем размер для определенного процессора... перебором


            int k = neighbourID[i] / S;
            int j = neighbourID[i] - k * S;
            int lenRecv = Compose_Q[k] * Compose_S[j];
            //printf("k=%i j=%i Recv len from %i: %i x %i = %i\n", k, j, neighbourID[i], Compose_Q[k], Compose_S[j], lenRecv);

            ///////////////////////////////////////////////////////////////////////
            AllRowRecv[i] = new double[lenRecv];

            if (this->id == neighbourID[i]) {
                double* Row = new double[lenRecv];
                int iter = 0;
                for (int j = 0; j < ScaleTable.wide; j++) {
                    for (int i = 0; i < ScaleTable.height; i++) {
                        Row[iter] = ScaleTable.Matrix[i][j];
                        iter++;
                    }
                }

                AllRowRecv[i] = Row;

                continue;
            }
                
            //прием
            AllRowRecv[i] = recv(lenRecv, neighbourID[i], 0);

        }
    }
    
    void multiply(int numproc, int* Compose_Q, int* Compose_S,  int Q, int S) {
        
        printf("\nУмножение\n");
        int columns = 0;
        for (int i = 0; i < Q; i++) {
            columns += Compose_Q[i];
        }

        int len = MainTable.height * columns;
        double* result = new double[len];

        //перемножаем и запоминаем суммы
        int r = 0;
        for (int i = 0; i < MainTable.height; i++) {

            for (int k = 0; k < Q; k++) { //цикл с столбцами умножения

                int len1 = MainTable.wide; //это можно как то заменить имея композиции и ID
                int delta = Compose_Q[k];

                for (int d = 0; d < delta; d++) {
                    result[r] = 0;

                    for (int j = 0; j < MainTable.wide; j++) {
                        double mult = MainTable.Matrix[i][j] * AllRowRecv[k][j + d * len1];
                        result[r] += mult;
                        //printf("\n(%i %i = %i): %f * %f = %f  result = %f", d, k, r, MainTable.Matrix[i][j], AllRowRecv[k][j+ d * len1], mult, result[r]);
                    }
                    printf("\n(%i) result = %f", r, result[r]);
                    r++;
                }
            }
            printf("\n");

        }
    
        //выведем
        printf("\nResult Array: ");
        for (int i = 0; i < MainTable.height * columns; i++) {
            printf("%f ", result[i]);
        }
        printf("\n");

        //отсылаем в последний для ссуммирования
        int cel = this->id / S;
        int sendToID = cel * S + (S - 1);

        send(result, len, sendToID, 1);

        if (this->id == sendToID) {
            double** resultRecv = new double* [S];
            sum = new double[len];

            //получаем
            for (int i = 0; i < S; i++) {

                int idFrom = this->id - (S - 1 - i);

                //printf("\n id=%i %i * %i = %i\n", idFrom, Compose_Q[k], columns, lenRecv);
                resultRecv[i] = new double[len];

                if (idFrom == this->id) resultRecv[i] = result;
                else resultRecv[i] = recv(len, idFrom, 1);  
            }
        
            //выводим
            printf("\nSum len=%i: ", len);
            for (int j = 0; j < len; j++) {
                sum[j] = 0;
                for (int i = 0; i < S; i++) {
                    sum[j] += resultRecv[i][j];
                }
                printf("%f ", sum[j]);
            }
            printf("\n");
            
            //отправим результат в самый последний процессор 
            send(sum, len, numproc-1, 2);
        }
    }
    
    void lustRecvResult(int numproc, int* Compose_Q, int Q, int S) {

        int columns = 0;
        for (int i = 0; i < Q; i++) {
            columns += Compose_Q[i];
        }

        //double** Result = new double* [Q];
        double** Result = new double* [columns];

        int offset = 0;
        //for (int i = 0; i < Q - 1; i++) {
        for (int i = 0; i < Q; i++) {
            int id = i * S + (S - 1);

            int k = id / S;
            int len = Compose_Q[k] * columns;

            double* recieve = new double[len];

            if (i < Q - 1) {
                recieve = recv(len, id, 2);
            }
            else {
                recieve = sum;
            }
           
            //остается это переписать в таблицу результатов и учесть Sum из функции multyply
            //офсет нужен для смещение при записи если в полученных результатах из других процессоров больше чем одна строках (так как передаем строкой всегда)
            for (int it = 0; it < Compose_Q[i]; it++) {
                //printf();
                Result[offset + it] = new double[columns];

                for (int j = 0; j < columns; j++) {
                    Result[offset + it][j] = recieve[j + columns*it];
                }
            }

            offset += Compose_Q[i];
        }
        printf("\nRESULT\n");
        for (int i = 0; i < columns; i++) {
            for (int j = 0; j < columns; j++) {
                printf("%f \t", Result[i][j]);
            }
            printf("\n");
        }

    }

    void send(double* RowSend, int len, int toID, int tag) {
        if (this->id == toID) return;
        printf("[%i]: Send to %i len = %i: ", this->id, toID, len);
        MPI_Send(RowSend, len, MPI_DOUBLE, toID, tag, MPI_COMM_WORLD);

        for (int i = 0; i < len; i++) printf("%f ", RowSend[i]); 
        printf("\n");
    }
    
    double* recv(int len, int fromID, int tag) {
        if (this->id == fromID) {
            printf("RECV this->id == fromID");
            system("pause");
        }
        double* Row = new double[len];
        printf("[%i]: Recv from %i len=%i: ", this->id, fromID, len);
        MPI_Recv(Row, len, MPI_DOUBLE, fromID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int k = 0; k < len; k++) {
            printf("%f\t ", Row[k]);
        }
        printf("\n");

        return Row;
    }

    double* recv1(int len, int fromID, int tag) {
        if (this->id == fromID) {
            printf("RECV this->id == fromID");
            system("pause");
        }
        double* Row = new double[len];
        printf("[%i]: Recv from %i len=%i: ", this->id, fromID, len);
        //MPI_Recv(Row, len, MPI_DOUBLE, fromID, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //for (int k = 0; k < len; k++) {
        //    printf("%f\t ", Row[k]);
        //}
        printf("\n");

        return Row;
    }
    
    void determine() {
        MainTable.clear();
        ScaleTable.clear();
    }
};

int main(int argc, char** argv){

	int numproc, id;
    int wide = 5, height = 5;
    int maxProc = wide * height - 1;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

 
    printf("\n====================================================\n");
    printf("\t\t\tPROC %i", id);
    printf("\n====================================================\n");
    //Конфиг должен быть в нулевом? там происходит выделения памяти для размеров таблиц... либо у каждого свой конфиг
	Config Config;
	Config.set(numproc, id, wide, height); //numproc, id, wide, hide
    int lustProc = Config.numproc - 1;
    if (id == 0) {
        //можно избавиться и заполнить у каждого процессора имея промежутки разбиения   
	    Table mainTable;
	    Table scaleTable;

	     mainTable.create(Config.wide, Config.height, Config.Compose_Q, Config.Compose_S);
	    scaleTable.create(Config.height, Config.wide, Config.Compose_S, Config.Compose_Q);

	     mainTable.printForMaintable();
	    scaleTable.printForMaintable();

	     mainTable.clear();
	    scaleTable.clear();
    }
    
    if (id < Config.numproc) {

        printf("==================Lust Proc = %i==================\n", id, lustProc);
        Proc Proces;
        //почти везде нужны Config.Compose_Q, Config.Compose_S, Config.Q, Config.S можно запомнить при создании
        Proces.set(id, Config.wideProc, Config.heightProc);
        Proces.fillTable(Config.Compose_Q, Config.Compose_S, Config.Q, Config.S);

        Proces.MainTable.print();
        printf("||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\n");
        Proces.ScaleTable.print();

        //обмен столбцами для умножений 
    
        Proces.calculateNeighbourID(Config.Q, Config.S);
        Proces.sendScaleRow(Config.numproc, Config.Compose_Q, Config.Compose_S, Config.Q, Config.S);
        Proces.recvScaleRow(Config.numproc, Config.Compose_Q, Config.Compose_S, Config.Q, Config.S);
        //MPI_Barrier(MPI_COMM_WORLD);

        //перемножение и отправка в последний в строке для ссумирования
        Proces.multiply(Config.numproc, Config.Compose_Q, Config.Compose_S, Config.Q, Config.S);


        if(id == Config.numproc-1) Proces.lustRecvResult(Config.numproc, Config.Compose_Q, Config.Q, Config.S);
        //MPI_Barrier(MPI_COMM_WORLD);
        Proces.determine();
    }

	MPI_Finalize();
    printf("\n====================================================\n");
	printf("\t\t\tEND %i", id);
    printf("\n====================================================\n");
}


