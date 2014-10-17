#include <iostream>
#include <string>
#include <sdsl/suffix_arrays.hpp>
#include <zmq.hpp>

using namespace std;
using namespace sdsl;

typedef struct cmdargs {
    string file;
    string port;
} cmdargs_t;

void
print_usage(string program)
{
    cout << program << " -f <file> -p <port>" << endl;
    cout << "where" << endl;
    cout << " -f <file>: the file which should be indexed and queried." << endl;
    cout << " -p <port>  : the port the daemon is running on." << endl;
};

cmdargs_t
parse_args(int argc,char* const argv[])
{
    cmdargs_t args;
    int op;
    args.file = "";
    args.port = "12345";
    while ((op=getopt(argc,argv,"f:p")) != -1) {
        switch (op) {
            case 'f':
                args.file = optarg;
                break;
            case 'p':
                args.port = optarg;
                break;
            default:
                print_usage(argv[0]);
        }
    }
    if (args.file == "") {
        cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

typedef csa_wt<wt_huff<rrr_vector<63>>,
        64, 64, text_order_sa_sampling<>,
        text_order_isa_sampling_support<>> csa_type;

int main(int argc, char* argv[])
{
    cmdargs_t args = parse_args(argc, argv);
    csa_type csa;
    string csa_file = args.file+".csa";
    if (!load_from_file(csa, csa_file)) {
        cout << "Index '" << csa_file << "' does not exist. Start construction." << endl;
        construct(csa, args.file, 1);
        cout << "Construction finished." << endl;
        store_to_file(csa, csa_file);
    }
    {
        cout << "Start daemon on port " << args.port << endl;
        zmq::context_t context(1);
        zmq::socket_t socket(context, ZMQ_REP);
        socket.bind(string("tcp://*:"+args.port).c_str());

        while (true) {
            zmq::message_t request;
            socket.recv(&request);
            char* req = (char*) request.data();
            uint64_t cnt = 0;
            string s(req+1, req+1+(uint8_t)req[0]);
            cnt = count(csa, s.begin(), s.end());
            zmq::message_t reply(sizeof(cnt));  // create a 8-byte message
            memcpy((void*) reply.data(), &cnt, sizeof(cnt));
            socket.send(reply);
        }

    }
}
