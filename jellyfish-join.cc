#include <string>
#include <fstream>
#include <iterator>

#include <boost/timer/timer.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <jellyfish/err.hpp>
#include <jellyfish/hash_counter.hpp>
#include <jellyfish/mer_iterator.hpp>
#include <jellyfish/mer_dna.hpp>
#include <jellyfish/jellyfish.hpp>
#include <jellyfish/large_hash_iterator.hpp>

using jellyfish::mer_dna;
using namespace std;
namespace po = boost::program_options;

struct fastq {
  string id;
  string seq;
  string plus;
  string qual;

  friend istream & operator>>(istream &is, fastq &record) {
    getline(is, record.id);
    getline(is, record.seq);
    getline(is, record.plus);
    getline(is, record.qual);
    return is;
  };
};

struct fasta {
  string id;
  string seq;

  friend istream & operator>>(istream &is, fasta &record) {
    getline(is, record.id);
    getline(is, record.seq);
    return is;
  };

};

template<typename Record, bool Canonical>
void join(istream_iterator<Record> iter, const mer_array *ary, const unsigned int kmer_length) {
  cerr << "=== Reading sequences ===" << endl;
  boost::timer::auto_cpu_timer t(cerr, 2);

  uint64_t val;
  mer_dna mer;
  Record record;
  istream_iterator<Record> EOS;
  unsigned int N;

  while(iter != EOS) {
    record = *iter++;
    N = record.seq.size() - kmer_length + 1;
    cout << record.id << endl;
    for(unsigned int i = 0; i < N; ++i) {
      // have to deal with 'N' values here (jellyfish deals with lowercases)
      if(!(mer.from_chars(record.seq.substr(i, kmer_length).c_str()))) {
        continue;
      }
      if(!ary->get_val_for_key(Canonical ? mer.get_canonical() : mer, &val)) {
        val = 0;
      }
      cout << val;
      if(i != N - 1) cout << "\t";
    }
    cout << endl;
  }
}

/* ========== main ========== */

int main(int argc, char *argv[]) {

  /* arguments */
  string jf_file;
  string fasta_file;
  string fastq_file;
  bool gzipped = false;

  po::options_description desc("take in a .jf matrix, and fasta/fastq fastq file, output for each kmer of each contig/read, its count in the jf file");
  desc.add_options()
    ("help,h", "Show help message")
    ("jf,j", po::value<string>(&jf_file)->required(), ".jf file")
    ("gzip,z", po::bool_switch(&gzipped), "fasta/fastq is gzipped")
    ("fasta", po::value<string>(&fasta_file), "fasta file")
    ("fastq", po::value<string>(&fastq_file), "fastq file");

  /* ---------- parse arguments ---------------- */
  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if(vm.count("help")) {
      cout << desc;
      exit(0);
    } 
    po::notify(vm);
  } catch(po::error &e) {
    cerr << "ERROR:" << e.what() << endl << desc;
    exit(1);
  }

  if(fasta_file == "" && fastq_file == "") {
    cerr << "ERROR: require one of fasta or fastq" << endl;
    exit(1);
  }

  /* ---------- Get header from jf file ----------- */
  ifstream jf_file_stream(jf_file);
  if(!jf_file_stream.good()) {
    cerr << "Error opening " << jf_file << endl;
    exit(1);
  }
  jellyfish::file_header header;
  header.read(jf_file_stream);

  uint64_t hash_size = header.size();
  unsigned int kmer_length = header.key_len() / 2;
  cerr << "hash_size: "   << hash_size   << endl 
       << "kmer_length: " << kmer_length << endl
       << "canonical?: "  << header.canonical() << endl;

  mer_dna::k(kmer_length);
  const int num_threads = 1;
  mer_hash hash(hash_size, header.key_len(), header.val_len(), num_threads, header.max_reprobe());
  hash.ary()->matrix(header.matrix());

  /* ---------- load the hash table with the k-mers from given file ---------- */
  {
    boost::timer::auto_cpu_timer t(cerr, 2);
    cerr << "=== Loading hash ===" << endl;
    binary_reader reader(jf_file_stream, &header);
    
    while(reader.next()) {
      hash.add(reader.key(), reader.val());
    }
    
  }

  bool is_fasta = fasta_file != "";
  ifstream ifs(is_fasta ? fasta_file : fastq_file, gzipped ? ios_base::binary : ios_base::in);
  boost::iostreams::filtering_istream fifs;
  if (gzipped) {
    fifs.push(boost::iostreams::gzip_decompressor());
  }
  fifs.push(ifs);

  if(is_fasta) {
    istream_iterator<fasta> iter(fifs);
    if(header.canonical())
      join<fasta, true>(iter, hash.ary(), kmer_length);
    else 
      join<fasta, false>(iter, hash.ary(), kmer_length);

  } else {
    istream_iterator<fastq> iter(fifs);
    if(header.canonical())
      join<fastq, true>(iter, hash.ary(), kmer_length);
    else 
      join<fastq, false>(iter, hash.ary(), kmer_length);
  }
 
  return 0;
}
