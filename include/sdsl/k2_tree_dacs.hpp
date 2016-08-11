/* Wrapper class around the dac c library */
/*
 * Adapted by Jan Broß
 */

#ifndef INCLUDE_COMPRESSION_DACS_H_
#define INCLUDE_COMPRESSION_DACS_H_

#include "../../external/dacs/include/dacs.h"
#include "../../external/dacs/include/bitrankw32int.h"
#include "k2_tree_helper.hpp"
#include "../../external/dacs/include/directcodes.h"


namespace sdsl {

    using namespace k2_treap_ns;

    class DAC {

    private:
        std::shared_ptr<sFTRep> m_comp_leaves;
    public:
        DAC() = default;

        struct Deleter {
            void operator()(sFTRep* p) const {
                ::destroyFT(p);
            }
        };

        /**
         * @param cnt Number of words in the vocabulary
         * @param size Size of each word in bytes
         */
        DAC(uint *list, uint listLength) : m_comp_leaves(::createFT(list, listLength),Deleter()) {
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            m_comp_leaves = std::shared_ptr<sFTRep>((sFTRep *) malloc(sizeof(struct sFTRep)),Deleter());

            read_member(m_comp_leaves->listLength, in);
            read_member(m_comp_leaves->nLevels, in);
            read_member(m_comp_leaves->tamCode, in);

            read_member(m_comp_leaves->tamtablebase, in);
            read_member(&m_comp_leaves->tablebase, m_comp_leaves->tamtablebase, in);
            //m_comp_leaves->base_bits = (ushort*) malloc(sizeof(ushort)*m_comp_leaves->nLevels);
            read_member(&m_comp_leaves->base_bits, m_comp_leaves->nLevels, in);
            read_member(&m_comp_leaves->base, m_comp_leaves->nLevels, in);
            read_member(&m_comp_leaves->levelsIndex, m_comp_leaves->nLevels + 1, in);
            read_member(&m_comp_leaves->iniLevel, m_comp_leaves->nLevels, in);
            read_member(&m_comp_leaves->rankLevels, m_comp_leaves->nLevels, in);
            read_member(&m_comp_leaves->levels, m_comp_leaves->tamCode / W + 1, in);

            load_bitrank(in);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += write_member(m_comp_leaves->listLength, out, child, "listLength");
            written_bytes += write_member(m_comp_leaves->nLevels, out, child, "nlevels");
            written_bytes += write_member(m_comp_leaves->tamCode, out, child, "tamcode");
            written_bytes += write_member(m_comp_leaves->tamtablebase, out, child, "tamtablebase");
            written_bytes += write_member(m_comp_leaves->tablebase, m_comp_leaves->tamtablebase, out, child, "tablebase");
            written_bytes += write_member(m_comp_leaves->base_bits, m_comp_leaves->nLevels, out, child, "basebits");
            written_bytes += write_member(m_comp_leaves->base, m_comp_leaves->nLevels, out, child, "base");
            written_bytes += write_member(m_comp_leaves->levelsIndex, m_comp_leaves->nLevels + 1, out, child,
                                          "levelsIdx");
            written_bytes += write_member(m_comp_leaves->iniLevel, m_comp_leaves->nLevels, out, child, "iniLevel");
            written_bytes += write_member(m_comp_leaves->rankLevels, m_comp_leaves->nLevels, out, child, "rankLevel");
            written_bytes += write_member(m_comp_leaves->levels, m_comp_leaves->tamCode / W + 1, out, child, "levels");

            written_bytes += serialize_bitrank(out, child);

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize_bitrank(std::ostream &out, structure_tree_node *v = nullptr) const {
            size_type written_bytes = 0;

            written_bytes += write_member(m_comp_leaves->bS->n, out, v, "bS_n");
            written_bytes += write_member(m_comp_leaves->bS->factor, out, v, "bS_factor");
            written_bytes += write_member(m_comp_leaves->bS->data, m_comp_leaves->bS->n / W + 1, out, v, "bs_data");
            written_bytes += write_member(m_comp_leaves->bS->Rs, m_comp_leaves->bS->n / m_comp_leaves->bS->s + 1, out,
                                          v, "bS_Rs");

            return written_bytes;

        }

        void load_bitrank(std::istream &in) {
            bitRankW32Int *bitRank = (bitRankW32Int *) malloc(sizeof(struct sbitRankW32Int));
            bitRank->b = 32;
            read_member(bitRank->n, in);
            bitRank->integers = bitRank->n/W;
            read_member(bitRank->factor, in);
            bitRank->s = bitRank->b * bitRank->factor;
            bitRank->owner = 1;
            read_member(&bitRank->data, bitRank->n / W + 1, in);
            read_member(&bitRank->Rs, bitRank->n / bitRank->s + 1, in);
            m_comp_leaves->bS = bitRank;
        }

        //! Swap operator
        void swap(DAC &dac) {
            if (this != &dac) {
                std::swap(m_comp_leaves, dac.m_comp_leaves);
            }
        }

        bool operator==(const DAC &rhs) const {
            if ((m_comp_leaves.get() == nullptr) && rhs.m_comp_leaves.get() == nullptr){
                    return true;
            }

            if ((m_comp_leaves.get() == nullptr) || rhs.m_comp_leaves.get() == nullptr){
                return false;
            }

            return (equalsFT(m_comp_leaves.get(), rhs.m_comp_leaves.get()));
        }

        void operator=(const DAC &dac) {
            m_comp_leaves = dac.m_comp_leaves;
        }

        uint accessFT(uint64_t i) const {
            //FIXME: no 64bit support
            return ::accessFT(m_comp_leaves.get(), (uint) i);
        }
    };

}  // namespace sdsl
#endif  // INCLUDE_COMPRESSION_VOCABULARY_H_
