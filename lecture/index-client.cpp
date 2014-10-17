#include <iostream>
#include <string>
#include <zmq.hpp>
#include <unistd.h>

using namespace std;

typedef struct cmdargs {
    std::string host;
} cmdargs_t;

void
print_usage(string program)
{
    cout << program << " -h <host>" << endl;
    cout << "where" << endl;
    cout << " -h <host>  : host of the daemon. Default: 127.0.0.1:12345" << endl;
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.host = "127.0.0.1:12345";
    while ((op=getopt(argc,argv,"p")) != -1) {
        switch (op) {
            case 'h':
                args.host = optarg;
                break;
            default:
                print_usage(argv[0]);
        }
    }
    return args;
}

struct query_t {
    uint64_t len;
    char* data;
};

char pattern[sizeof(uint8_t)+1024];

int main(int argc, char* argv[])
{
    cmdargs_t args = parse_args(argc, argv);

    zmq::context_t context(1);
    zmq::socket_t socket(context, ZMQ_REQ);
    socket.connect(std::string("tcp://"+args.host).c_str());
    if (!socket.connected()) {
        cerr << "Connecting to daemon failed." << endl;
        return 1;
    }


    char* pat = pattern+8;


    while (cin.getline(pat, 1024)) {
        uint8_t len = strlen(pat);
        zmq::message_t request(sizeof(len)+len);
        memcpy((void*)request.data(), &len, sizeof(len));
        char* rp = ((char*)request.data())+sizeof(len);
        memcpy((void*)rp, pat, len);
        socket.send(request);
        zmq::message_t output;
        socket.recv(&output);
        uint64_t cnt=0;
        memcpy((void*)&cnt, output.data(), sizeof(cnt));
        cout<<cnt<<endl;
    }
}
