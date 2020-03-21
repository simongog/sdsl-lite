#include <sdsl/int_vector.hpp>
#include <sdsl/int_vector_mapper.hpp>
#include <iostream>
#include <iomanip>

#include "utils.hpp"
#include "constants.hpp"
#include "collection.hpp"
#include "logging.hpp"

typedef struct cmdargs {
    std::string input_file;
    uint8_t     num_bytes;
    std::string collection_dir;
} cmdargs_t;

void print_usage(const char* program)
{
    fprintf(stdout, "%s -i <input file> [-n <input num bytes>] -c <col dir>\n", program);
    fprintf(stdout, "where\n");
    fprintf(stdout, "  -i <input file>      : the input file.\n");
    fprintf(stdout, "  -n <input num bytes> : #bytes per symbol or 0 for int_vector<0>. (default: 1)\n");
    fprintf(stdout, "  -c <collection dir>  : the collection dir.\n");
};

cmdargs_t parse_args(int argc, const char* argv[])
{
    cmdargs_t args;
    int op;
    args.input_file = "";
    args.collection_dir = "";
    args.num_bytes = 1;
    while ((op = getopt(argc, (char* const*)argv, "i:n:c:")) != -1) {
        switch (op) {
            case 'i':
                args.input_file = optarg;
                break;
            case 'n':
                args.num_bytes = std::stoi(optarg);
                break;
            case 'c':
                args.collection_dir = optarg;
                break;
        }
    }
    if (args.collection_dir == "" || args.input_file == "") {
        std::cerr << "Missing command line parameters.\n";
        print_usage(argv[0]);
        exit(EXIT_FAILURE);
    }
    return args;
}

int main(int argc, const char* argv[])
{
    log::start_log(argc, argv);
    cmdargs_t args = parse_args(argc, argv);


    /* (1) create collection dir */
    LOG(INFO) << "Creating collection directory structure in '" << args.collection_dir << "'";
    utils::create_directory(args.collection_dir);

    /* (2) create sub dirs */
    utils::create_directory(args.collection_dir+"/tmp");
    utils::create_directory(args.collection_dir+"/index");
    utils::create_directory(args.collection_dir+"/patterns");
    utils::create_directory(args.collection_dir+"/results");

    if (args.num_bytes == 1) {
        /* (3) check if input file exists */
        std::ifstream ifs(args.input_file,std::ios::binary);
        if (!ifs) {
            LOG(FATAL) << "Error opening input file '" << args.input_file << "'";
        }
        ifs >> std::noskipws; // dont skip white spaces!

        /* (4) copy file to sdsl format */
        auto buf = sdsl::write_out_buffer<0>::create(args.collection_dir+"/"+ consts::KEY_PREFIX + consts::KEY_TEXT);
        {
            gm_timer tm("COPY TO SDSL FORMAT",true);
            std::copy(std::istream_iterator<uint8_t>(ifs),
                      std::istream_iterator<uint8_t>(),
                      std::back_inserter(buf));
        }
        {
            for (size_t i=0; i<buf.size(); ++i) {
                if (buf[i] == 0 or buf[i] == '\n' or buf[i] == '\r' or buf[i] == '\f')
                    buf[i] = ' ';
            }
            gm_timer tm("BIT COMPRESS",true);
            sdsl::util::bit_compress(buf);
        }

        LOG(INFO) << "num_syms = " << buf.size();
        LOG(INFO) << "log2(sigma) = " << (int) buf.width();
        LOG(INFO) << "text size = " << (buf.width()*buf.size())/(8*1024*1024) << " MiB";
    } else { // handle integer alphabets
        sdsl::int_vector<> buf;
        if (!sdsl::load_vector_from_file(buf, args.input_file, args.num_bytes)) {
            LOG(FATAL) << "Error opening input file '" << args.input_file << "'";
        }
        for (size_t i=0; i<buf.size(); ++i)
            if (buf[i] == 0)
                buf.resize(i);
        sdsl::util::bit_compress(buf);
        sdsl::store_to_file(buf, args.collection_dir+"/"+ consts::KEY_PREFIX + consts::KEY_TEXT);

        LOG(INFO) << "num_syms = " << buf.size();
        LOG(INFO) << "log2(sigma) = " << (int) buf.width();
        LOG(INFO) << "text size = " << (buf.width()*buf.size())/(8*1024*1024) << " MiB";
    }

    return EXIT_SUCCESS;
}
