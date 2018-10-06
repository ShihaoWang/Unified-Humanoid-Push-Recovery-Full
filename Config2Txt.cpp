// g++ -Wall -O3 -o Config2Txt Config2Txt.cpp
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

int not_main()
{
    string dataArray[36];
    ifstream input_file;
    input_file.open("robot_angle_init.config");
    if(input_file.is_open())
    {
        char temp;
        input_file>>temp;
        input_file>>temp;
        for(int i = 0; i<36; i++)
        {
            input_file>>dataArray[i];
        }
    }
    input_file.close();

    // After the loading, this program will rewrite them into a column fashion
    ofstream output_file;
    int Act_Link_Ind_arr[] = {0, 2, 4, 8, 9, 10, 14, 15, 16, 22, 25, 29, 32};
    int n = sizeof(Act_Link_Ind_arr)/sizeof(Act_Link_Ind_arr[0]);

    std::vector<int> Act_Link_Ind(Act_Link_Ind_arr, Act_Link_Ind_arr + n);

    output_file.open("robot_angle_init.txt", std::ofstream::out);
    for (int i = 0; i < 13; i++)
    {
        output_file<<dataArray[Act_Link_Ind[i]]<<endl;
    }
    output_file.close();
    return 0;
}
