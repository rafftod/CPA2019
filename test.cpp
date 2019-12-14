#include <iostream>


#ifdef __linux__
#include <unistd.h>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __APPLE__
#include <unistd.h>
#endif

int check_cores()
{
    int numCPU = 0;
    #ifdef _WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    numCPU = sysinfo.dwNumberOfProcessors
    #endif

    #ifdef __linux__
    numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    #endif

    #ifdef __APPLE__
    numCPU = sysconf(_SC_NPROCESSORS_ONLN);
    #endif
}

int main(int argc, char const *argv[])
{
    int num = check_cores();
    if (num != 0)
    {
        std::cout<<"This machine has: "<<num<<" cores"<<std::endl;
    }
    
    return 0;
}
