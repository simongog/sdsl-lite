/*
 * ----------------------------------------------------------------------------
 * "THE BEER-WARE LICENSE" (Revision 42):
 * <nlehmann@dcc.uchile.cl> wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer in return Nicolás Lehmann
 * ----------------------------------------------------------------------------
 */

/* Wavelet Tree over Dense Code. --
   A word index using wavelet tree strategy over compressed text.

   Programmed by Susana Ladra.

   Author's contact: Susana Ladra, Databases Lab, University of
   A Coruna. Campus de Elvina s/n. Spain  sladra@udc.es

   Copyright (C) 2007  Nieves R. Brisaboa, Antonio Farina and Susana Ladra
   Author's contact: susanaladra@gmail.com

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
   */

/*
 * Adapted by Jan Broß
 */

#ifndef INCLUDE_COMPRESSION_DACS_H_
#define INCLUDE_COMPRESSION_DACS_H_

#include <algorithm>
#include <fstream>
#include "k2_tree_bitrankw32int.h"
#include "k2_tree_helper.hpp"
#include "k2_tree_dac_helper.hpp"



namespace sdsl {

    using namespace k2_treap_ns;

    class DAC {

    private:
        const ushort FACT_RANK = 20;

        uint m_listLength;
        byte m_nLevels;
        uint m_tamCode;
        uint *m_levels;
        uint *m_levelsIndex;
        uint *m_iniLevel;
        uint *m_rankLevels;
        bit_rank_w32int m_bS;
        //uint * bits_bitmap;
        uint *m_base;
        ushort *m_base_bits;
        uint *m_tablebase;
        uint m_tamtablebase;
    public:
        DAC() = default;

        /**
         * @param cnt Number of words in the vocabulary
         * @param size Size of each word in bytes
         */
        DAC(uint *list, uint listLength) : m_listLength(listLength) {
            uint *levelSizeAux;
            uint *cont;
            uint *contB;

            ushort *kvalues;
            uint nkvalues;

            register uint i;
            uint l;
            uint value, newvalue;
            uint bits_BS_len = 0;

            //ushort kvalues[4] = {0,2,4,8};
            //uint nkvalues=4;

            uint maxInt = 0;

            for (i = 0; i < listLength; i++)
                if (maxInt < list[i])
                    maxInt = list[i];

            //fprintf(stderr,"maxInt : %d\n",maxInt);

            maxInt++;

            uint *weight = (uint *) malloc(sizeof(uint) * maxInt);


            for (l = 0; l < maxInt; l++)
                weight[l] = 0;

            for (i = 0; i < listLength; i++)
                weight[list[i]]++;


            uint *acumFreq = (uint *) malloc(sizeof(uint) * (maxInt + 1));

            acumFreq[0] = 0;
            for (i = 0; i < maxInt; i++)
                acumFreq[i + 1] = acumFreq[i] + weight[i];

            free(weight);


            kvalues = optimizationk(acumFreq, maxInt, &nkvalues);


            free(acumFreq);


            ushort kval;
            uint oldval = 0;
            uint newval = 0;

            i = 0;
            uint multval = 1;
            do {
                oldval = newval;
                if (i >= nkvalues) {
                    kval = 1 << (kvalues[nkvalues - 1]);
                }
                else
                    kval = 1 << (kvalues[i]);
                multval *= kval;
                newval = oldval + multval;//mypow(kval,i);
                i++;
            }
            while (oldval < newval);

            m_tamtablebase = i;
            m_tablebase = (uint *) malloc(sizeof(uint) * m_tamtablebase);
            levelSizeAux = (uint *) malloc(sizeof(uint) * m_tamtablebase);
            cont = (uint *) malloc(sizeof(uint) * m_tamtablebase);
            contB = (uint *) malloc(sizeof(uint) * m_tamtablebase);

            oldval = 0;
            newval = 0;
            multval = 1;
            for (i = 0; i < m_tamtablebase; i++) {
                oldval = newval;
                if (i >= nkvalues) {
                    kval = 1 << (kvalues[nkvalues - 1]);
                }
                else
                    kval = 1 << (kvalues[i]);
                multval *= kval;
                newval = oldval + multval;//mypow(kval,i);
                m_tablebase[i] = oldval;
            }


            for (i = 0; i < m_tamtablebase; i++) {
                levelSizeAux[i] = 0;

            }
            //Reservando espacio para los niveles

            for (i = 0; i < listLength; i++) {
                value = list[i];
                for (uint j = 0; j < m_tamtablebase; j++)
                    if (value >= m_tablebase[j])
                        levelSizeAux[j]++;
            }
            //
            uint o = 0;
            //	//Contadores
            while ((o < m_tamtablebase) && (levelSizeAux[o] != 0)) {
                o++;
            }
            m_nLevels = o;
            m_levelsIndex = (uint *) malloc(sizeof(uint) * (m_nLevels + 1));
            bits_BS_len = 0;

            m_base = (uint *) malloc(sizeof(uint) * m_nLevels);
            m_base_bits = (ushort *) malloc(sizeof(ushort) * m_nLevels);

            for (i = 0; i < m_nLevels; i++) {
                if (i >= nkvalues) {
                    m_base[i] = 1 << (kvalues[nkvalues - 1]);
                    m_base_bits[i] = kvalues[nkvalues - 1];
                }
                else {
                    m_base[i] = 1 << (kvalues[i]);
                    m_base_bits[i] = kvalues[i];
                }
            }

            uint tamLevels = 0;


            tamLevels = 0;
            for (i = 0; i < m_nLevels; i++)
                tamLevels += m_base_bits[i] * levelSizeAux[i];

            m_iniLevel = (uint *) malloc(sizeof(uint) * m_nLevels);
            m_tamCode = tamLevels;
            uint indexLevel = 0;
            m_levelsIndex[0] = 0;
            for (uint j = 0; j < m_nLevels; j++) {
                m_levelsIndex[j + 1] = m_levelsIndex[j] + levelSizeAux[j];
                //if(j>0){
                m_iniLevel[j] = indexLevel;
                cont[j] = m_iniLevel[j];
                indexLevel += levelSizeAux[j] * m_base_bits[j];
                //}
                contB[j] = m_levelsIndex[j];

            }


            m_levels = (uint *) malloc(sizeof(uint) * (tamLevels / W + 1));
            //fprintf(stderr,"tamaño de levels: %d\n",sizeof(uint)*(tamLevels/W+1));
            bits_BS_len = m_levelsIndex[m_nLevels - 1] + 1;
            //Se pone el ultimo a 0 para ahorrarnos comparaciones despues en la descompresion
            uint *bits_BS = (uint *) malloc(sizeof(uint) * (bits_BS_len / W + 1));
            int j,k;
            for (i = 0; i < ((bits_BS_len) / W + 1); i++)
                bits_BS[i] = 0;
            for (i = 0; i < listLength; i++) {
                value = list[i];
                j = m_nLevels - 1;

                while (j >= 0) {
                    if (value >= m_tablebase[j]) {

                        newvalue = value - m_tablebase[j];

                        for (k = 0; k < j; k++) {


                            k2_dac_helper::bitwrite(m_levels, cont[k], m_base_bits[k], newvalue % m_base[k]);
                            cont[k] += m_base_bits[k];
                            contB[k]++;
                            newvalue = newvalue / m_base[k];
                        }
                        k = j;


                        k2_dac_helper::bitwrite(m_levels, cont[j], m_base_bits[j], newvalue % m_base[j]);
                        cont[j] += m_base_bits[j];
                        contB[j]++;
                        if (j < m_nLevels - 1) {
                            bitset(bits_BS, contB[j] - 1);
                        }

                        break;
                    }
                    j--;
                }
                //Para j=0 solo se cubre el bitmap

            }
            //Para simular ultimo array:

            bitset(bits_BS, bits_BS_len - 1);
            //m_bits_bitmap = bits_BS;
            m_bS = bit_rank_w32int(bits_BS, bits_BS_len, 1, 20);

            m_rankLevels = (uint *) malloc(sizeof(uint) * m_nLevels);
            for (i = 0; i < m_nLevels; i++)
                m_rankLevels[i] = m_bS.rank(m_levelsIndex[i] - 1);


            free(cont);
            free(contB);
            free(levelSizeAux);
            free(kvalues);
        }

        //! Loads the data structure from the given istream.
        void load(std::istream &in) {
            read_member(m_listLength, in);
            read_member(m_nLevels, in);
            read_member(m_tamCode, in);

            read_member(m_tamtablebase, in);
            read_member(m_base_bits, m_nLevels, in);
            read_member(m_base, m_nLevels, in);
            read_member(m_levelsIndex, m_nLevels+1, in);
            read_member(m_iniLevel, m_nLevels, in);
            read_member(m_rankLevels, m_nLevels, in);
            read_member(m_levels, m_tamCode/W+1, in);

            m_bS.load(in);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
                            std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;

            written_bytes += write_member(m_listLength, out, child, "listLength");
            written_bytes += write_member(m_nLevels, out, child, "nlevels");
            written_bytes += write_member(m_tamCode, out, child, "tamcode");

            written_bytes += write_member(m_tamtablebase, out, child, "tamtablebase");

            written_bytes += write_member(m_base_bits, m_nLevels, out, child, "basebits");
            written_bytes += write_member(m_base, m_nLevels, out, child, "base");
            written_bytes += write_member(m_levelsIndex, m_nLevels+1, out, child, "levelsIdx");
            written_bytes += write_member(m_iniLevel, m_nLevels, out, child, "iniLevel");
            written_bytes += write_member(m_rankLevels, m_nLevels, out, child, "rankLevel");
            written_bytes += write_member(m_levels, m_tamCode/W+1, out, child, "levels");

            written_bytes += m_bS.serialize(out, child, "bitrank32W");

            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        ~DAC() {
            free(m_levelsIndex);
            free(m_iniLevel);
            free(m_rankLevels);
            free(m_tablebase);
            free(m_levels);
            //destroyBitmap(m_bS);
            //free(m_bits_bitmap);
            free(m_base);
            free(m_base_bits);
        }

        //! Swap operator
        void swap(DAC &dac) {
            if (this != &dac) {
                std::swap(m_listLength, dac.m_listLength);
                std::swap(m_nLevels, dac.m_nLevels);
                std::swap(m_tamCode, dac.m_tamCode);
                std::swap(m_tamtablebase, dac.m_tamtablebase);
                std::swap(m_base_bits, dac.m_base_bits);
                std::swap(m_base, dac.m_base);
                std::swap(m_levelsIndex, dac.m_levelsIndex);
                std::swap(m_iniLevel, dac.m_iniLevel);
                std::swap(m_rankLevels, dac.m_rankLevels);
                std::swap(m_levels, dac.m_levels);
            }
        }

/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */

        uint accessFT(uint param) const {
            uint mult = 0;
            //register uint i;
            register uint j;
            //uint * rankLevel = m_rankLevels;
            //byte * list;
            //uint n , partialSum=0, sumAux=0;
            uint partialSum = 0;
            //uint ini = param-1;
            uint ini = param;
            //bitRankW32Int * bS = m_bS;
            //uint * bsData = m_bS->data;
            uint nLevels = m_nLevels;
            //uint levelIndex;
            uint *level;
            uint readByte;
            uint cont, pos, rankini;


            //	fprintf(stderr,"Queriendo leer la posicion: %d\n",ini);
            partialSum = 0;
            j = 0;
            level = m_levels;//+ (m_levelsIndex[j]>>1);
            //		fprintf(stderr,"primera posicion del array levels: %d %d\n",level[0],(byte)level[0]);
            //cont=ini+(m_levelsIndex[j]&0x1);
            pos = m_levelsIndex[j] + ini;

            mult = 0;
            //readByte = ((*(level+(cont>>1)))>>(BASE_BITS*(cont&0x1)))&0xF;
            cont = m_iniLevel[j] + ini * m_base_bits[j];
            //		fprintf(stderr,"leyendo los %d bits \n",m_base_bits[j]);
            //		fprintf(stderr,"inilevel: %d \n",m_iniLevel[j]);
            //		fprintf(stderr,"ini: %d\n",ini);
            //		fprintf(stderr,"base_bits: %d\n",m_base_bits[j]);
            //
            //		fprintf(stderr,"%d posicion\n",cont);
            //		fprintf(stderr,"nivel %d\n",j);
            //		fprintf(stderr,"leyendo los %d bits desde %d posicion del nivel %d\n",m_base_bits[j],cont,j);
            readByte = k2_dac_helper::bitread(level, cont, m_base_bits[j]);
            //		fprintf(stderr,"readbyte: %d\n",readByte);
            //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);
            //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level[j]+(cont[j]>>1))));
            //	fprintf(stderr,"readByte: %d\n",readByte);
            //fprintf(stderr,"pos[%d]= %d\n",j,pos[j]);
            if (nLevels == 1) {
                return readByte;
            }
            while ((!bitget(m_bS.data, pos))) {
                //fprintf(stderr,"pos[%d]= %d\n",j,pos);
                rankini = m_bS.rank(m_levelsIndex[j] + ini - 1) - m_rankLevels[j];
                ini = ini - rankini;

                //fprintf(stderr,"readByte: %d\n",readByte);
                //partialSum = partialSum+ readByte*pot256[j];
                partialSum = partialSum + (readByte << mult);

                mult += m_base_bits[j];
                j++;
                //fprintf(stderr,"pos[%d] = %d lsit-> levels %d cont[%d] = %d \n",j-1,pos[j-1]-1,m_levels, j,cont[j]);

                //level=m_levels + (m_levelsIndex[j]>>1);
                //cont=ini+(m_levelsIndex[j]&0x1);
                cont = m_iniLevel[j] + ini * m_base_bits[j];
                pos = m_levelsIndex[j] + ini;

                //readByte = ((*(level+(cont>>1)))>>(BASE_BITS*(cont&0x1)))&0xF;
                //		fprintf(stderr,"leyendo los %d bits desde %d posicion del nivel %d\n",m_base_bits[j],cont,j);
                readByte = k2_dac_helper::bitread(level, cont, m_base_bits[j]);
                //		fprintf(stderr,"readbyte: %d\n",readByte);
                //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);

                //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level+(cont>>1))));
                if (j == nLevels - 1) {
                    break;
                }


            }
            //fprintf(stderr,"readByte: %d\n",readByte);
            //partialSum = partialSum + readByte*pot256[j] + m_tablebase[j];
            partialSum = partialSum + (readByte << mult) + m_tablebase[j];
            //fprintf(stderr,"partialSum = %d, param = %u \n",partialSum,param);
            //list[i]=sum
            //	fprintf(stderr,"devolviendo %d\n",partialSum);
            return partialSum;

        }


/*-----------------------------------------------------------------

  ---------------------------------------------------------------- */

        uint *decompressFT(uint n) {
            uint mult = 0;
            register uint i;
            register uint j;
            //uint partialSum=0, sumAux=0;
            uint partialSum = 0;
            //uint ini = 0;
            //bitRankW32Int * bS = m_bS;
            //uint * bsData = m_bS->data;
            uint nLevels = m_nLevels;
            //uint levelIndex;
            uint *level = m_levels;
            byte readByte;
            uint *list = (uint *) malloc(sizeof(uint) * n);
            uint *cont = (uint *) malloc(sizeof(byte *) * m_nLevels);
            uint *pos = (uint *) malloc(sizeof(uint) * m_nLevels);

            for (j = 0; j < nLevels; j++) {
                cont[j] = m_iniLevel[j];
                pos[j] = m_levelsIndex[j];
                //fprintf(stderr,"cont[%d] = %d, pos[%d]= %d\n",j,cont[j],j,pos[j]);
            }
            // uint * rankLevel = m_rankLevels;
            // byte * list;

            //uint cont,pos, rankini;
            //fprintf(stderr,"Queriendo leer la posicion: %d\n",ini);
            for (i = 0; i < n; i++) {
                partialSum = 0;
                j = 0;
                //level=m_levels + (m_levelsIndex[j]>>1);
                //cont=ini+(m_levelsIndex[j]&0x1);
                //pos=m_levelsIndex[j]+ini;

                mult = 0;
                readByte = k2_dac_helper::bitread(level, cont[j], m_base_bits[j]);
                cont[j] += m_base_bits[j];
                //readByte = ((*(level+(cont[j]>>1)))>>(BASE_BITS*(cont[j]&0x1)))&0xF;
                //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);
                //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level[j]+(cont[j]>>1))));

                //fprintf(stderr,"pos[%d]= %d\n",j,pos[j]);
                while ((!bitget(m_bS.data, pos[j]))) {
                    //fprintf(stderr,"pos[%d]= %d\n",j,pos[j]);
                    //rankini = rank(m_bS, m_levelsIndex[j]+ini-1) - m_rankLevels[j];
                    //ini = ini-rankini;
                    pos[j]++;
                    //cont[j]++;
                    //fprintf(stderr,"readByte: %d\n",readByte);
                    //partialSum = partialSum+ readByte*pot256[j];
                    partialSum = partialSum + (readByte << mult);
                    mult += m_base_bits[j];
                    j++;

                    //fprintf(stderr,"pos[%d] = %d lsit-> levels %d cont[%d] = %d \n",j-1,pos[j-1]-1,m_levels, j,cont[j]);


                    //level=m_levels + (m_levelsIndex[j]>>1);
                    //cont=ini+(m_levelsIndex[j]&0x1);
                    //pos=m_levelsIndex[j]+ini;
                    readByte = k2_dac_helper::bitread(level, cont[j], m_base_bits[j]);
                    cont[j] += m_base_bits[j];
                    //readByte = ((*(level+(cont[j]>>1)))>>(BASE_BITS*(cont[j]&0x1)))&0xF;
                    //fprintf(stderr,"contenido del vector: %d readByte: %d con cont[%d]: %d\n",*(level[j]+(cont[j]>>1)),readByte,j,cont[j]);

                    //fprintf(stderr,"readByte: %d... %d\n",readByte,(*(level[j]+(cont[j]>>1))));
                    if (j == nLevels - 1) {
                        break;
                    }

                }
                //fprintf(stderr,"readByte: %d\n",readByte);
                //partialSum = partialSum + readByte*pot256[j] + m_tablebase[j];
                if (j < nLevels - 1) {
                    pos[j]++;
                }


                partialSum = partialSum + (readByte << mult) + m_tablebase[j];
                //fprintf(stderr,"sum = %d, partialSum = %d, param = %u \n",sum,partialSum,param);
                //list[i]=sum
                //if(i==16404)
                //fprintf(stderr,"devolviendo %d\n",partialSum);
                list[i] = partialSum;
            }
            free(cont);
            free(pos);
            return list;
        }

        uint getListLength() {
            return m_listLength;
        }

        bool operator==(const DAC &rhs) const{
            if (m_listLength != rhs.m_listLength || m_nLevels != rhs.m_nLevels ||
                    m_tamCode != rhs.m_tamCode || m_tamtablebase != rhs.m_tamtablebase)
                return false;

            for (uint i = 0; i < m_nLevels; ++i) {
                if (m_base_bits[i] != rhs.m_base_bits[i] ||
                    m_base[i] != rhs.m_base[i] ||
                    m_iniLevel[i] != rhs.m_iniLevel[i] ||
                    m_rankLevels[i] != rhs.m_rankLevels[i] ||
                    m_levelsIndex[i] != rhs.m_levelsIndex[i])
                    return false;
            }

            if (m_levelsIndex[m_nLevels] != rhs.m_levelsIndex[m_nLevels])
                return false;

            for (uint i = 0; i < m_tamtablebase; ++i)
                if (m_tablebase[i] != rhs.m_tablebase[i])
                    return false;

            for (uint i = 0; i < m_tamCode / W + 1; ++i)
                if (m_levels[i] != rhs.m_levels[i])
                    return false;
            return (m_bS == rhs.m_bS);
        }

        void operator=(const DAC &dac){
            m_listLength = dac.m_listLength;
            m_nLevels = dac.m_nLevels;
            m_tamCode = dac.m_tamCode;
            m_tamtablebase = dac.m_tamtablebase;
            m_base_bits = dac.m_base_bits;
            m_base = dac.m_base;
            m_levelsIndex = dac.m_levelsIndex;
            m_iniLevel = dac.m_iniLevel;
            m_rankLevels = dac.m_rankLevels;
            m_levels = dac.m_levels;
        }

    private:
        uint mypow(uint base, uint exponente) {
            uint result = 1;
            for (uint i = 0; i < exponente; i++) {
                result *= base;
            }
            return result;
        }


        ushort *optimizationk(uint *acumFreqs, int maxInt, uint *nkvalues) {
            int sizeVoc = maxInt;

            //uint listLength = acumFreqs[sizeVoc];
            uint nBits = k2_dac_helper::bits(sizeVoc);

            //fprintf(stderr,"Length of the list: %d, max bits: %d\n",listLength,nBits);

            //Esta tabla tiene el tamanho que ocupa la mejor opcion para los primeros x bits
            long *tableSize = (long *) malloc(sizeof(long) * (nBits + 1));

            //Estos dos arrays indican los bs optimos hasta este punto
            //El primer array dice cuantos niveles hay hasta ese punto y el otro los bits en cada nivel
            uint *tableNLevels = (uint *) malloc(sizeof(uint) * (nBits + 1));
            uint **tableKvalues = (uint **) malloc(sizeof(uint *) * (nBits + 1));

            uint j;
            ulong maxSize = 0, maxPos = 0;
            int posVocInf, posVocSup;
            ulong currentSize;
            //Para la optimizacion, hay que mirar el máximo de todas las opciones anteriores anhadiendo
            //hasta el bit actual

            tableSize[0] = 0;
            tableNLevels[0] = 0;
            tableKvalues[0] = NULL;


            uint i;
            for (i = 1; i <= nBits; i++) {
                maxSize = -1;
                maxPos = 0;


                for (j = 0; j < i; j++) {
                    if (i == nBits)
                        posVocInf = 0;
                    else
                        posVocInf = 1 << (nBits - i);
                    posVocSup = (1 << (nBits - j));
                    if (posVocSup >= sizeVoc)
                        posVocSup = sizeVoc;
                    if (j == 0)
                        currentSize = tableSize[j] + ((ulong) (acumFreqs[sizeVoc] - acumFreqs[posVocInf])) * ((i - j));
                    else
                        currentSize =
                                tableSize[j] + ((ulong) (acumFreqs[sizeVoc] - acumFreqs[posVocInf])) * ((i - j) + 1) +
                                (acumFreqs[sizeVoc] - acumFreqs[posVocInf]) / FACT_RANK;

                    if (maxSize > currentSize) {
                        maxSize = currentSize;
                        maxPos = j;
                    }

                }

                tableSize[i] = maxSize;
                tableNLevels[i] = tableNLevels[maxPos] + 1;
                tableKvalues[i] = (uint *) malloc(sizeof(uint) * tableNLevels[i]);
                for (j = 0; j < tableNLevels[i] - 1; j++)
                    tableKvalues[i][j] = tableKvalues[maxPos][j];
                tableKvalues[i][tableNLevels[i] - 1] = i - maxPos;


            }


            ulong sumaTotal = 0;
            int bitCountInf = 0, bitCountSup = 0, bitsCount;
            for (j = 0; j < tableNLevels[nBits]; j++) {
                bitsCount = tableKvalues[nBits][tableNLevels[nBits] - 1 - j];
                bitCountSup += bitsCount;
                if (bitCountInf == 0)
                    posVocInf = 0;
                else
                    posVocInf = 1 << bitCountInf;
                posVocSup = (1 << bitCountSup);
                if (posVocSup >= sizeVoc)
                    posVocSup = sizeVoc;
                if (j == tableNLevels[nBits])
                    sumaTotal += ((ulong) (acumFreqs[sizeVoc] - acumFreqs[posVocInf])) * bitsCount;

                else
                    sumaTotal += ((ulong) (acumFreqs[sizeVoc] - acumFreqs[posVocInf])) * (bitsCount + 1) +
                                 (acumFreqs[sizeVoc] - acumFreqs[posVocInf]) / FACT_RANK;


                bitCountInf += bitsCount;
            }


            (*nkvalues) = tableNLevels[nBits];

            ushort *kvalues = (ushort *) malloc(sizeof(ushort) * tableNLevels[nBits]);

            for (j = 0; j < tableNLevels[nBits]; j++) {
                bitsCount = tableKvalues[nBits][tableNLevels[nBits] - 1 - j];

                kvalues[j] = bitsCount;


                bitCountInf += bitsCount;
            }


            free(tableSize);
            for (i = 1; i <= nBits; ++i)
                free(tableKvalues[i]);
            free(tableKvalues);
            free(tableNLevels);
            return kvalues;
        }
    };

}  // namespace sdsl
#endif  // INCLUDE_COMPRESSION_VOCABULARY_H_
