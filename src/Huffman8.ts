/*
 * Mostly taken from huffman.c: http://www.romhacking.net/utilities/826/
 */
import BitConverter from './BitConverter'

class HuffmanNode {
    symbol: number;
    weight: number;
    leafs: number;
    
    dad: HuffmanNode;
    lson: HuffmanNode;
    rson: HuffmanNode;
}

class HuffmanCode {
    nbits: number;
    codework: Array<number>;
}

export default class Huffman8 {
    const CMD_CODE: number = 0x28; // 8-bit Huffman magic number
    const HUF_LNODE: number = 0;
    const HUF_RNODE: number = 1;

    const HUF_SHIFT: number = 1;
    const HUF_MASK: number = 0x80;
    const HUF_MASK4: number = 0x80000000;

    const HUF_LCHAR: number = 0x80;
    const HUF_RCHAR: number = 0x40;
    const HUF_NEXT: number = 0x3F;

    const HUF_MAXIM: number = 0x1400000;

    const max_symbols: number = 0x100;

    tree: Array<HuffmanNode>;
    codetree: Array<number>;
    codemask: Array<number>;
    codes: Array<HuffmanCode>;
    num_leafs: number;
    num_nodes: number;
    freqs: Array<number>;
    num_bits: number;
    
    public decompress(data: Array<number>) {
        let header = BitConverter.toUInt32(data, 0);
        num_bits = header & 0xF;
        let uncompressedLength = header >> 8;
        let uncompressed = new Array<number>(uncompressedLength);
        let pak_pos = 4;
        pak_pos+= (data[pak_pos] + 1) << 1;
        let raw_pos = 0;
        const tree_ofs = 4;
        let nbits = 0;
        let pos = data[tree_ofs + 1];
        let next = 0;
        let mask4 = 0;
        let code = BitConverter.toUInt32(data, pak_pos);
        
        while (raw_pos < uncompressed.length) {
            if((mask >>= HUF_SHIFT) === 0) {
                if(pak_pos + 3 >= data.length) {
                    break;
                }
                code = BitConverter.toUInt32(data, pak_pos);
                pak_pos+= 4;
                mask4 = HUF_MASK4;
            }
            
            next+= ((pos & HUF_NEXT) + 1) << 1;
            
            let ch;
            if ((code & mask4) === 0) {
                ch = pos & HUF_LCHAR;
                pos = data[tree_ofs + next];
            } else {
                ch = pos & HUF_RCHAR;
                pos = data[tree_ofs + next + 1];
            }
            
            if (ch != 0)
            {
                uncompressed[raw_pos] |= (pos << nbits);
                nbits = (nbits + num_bits) & 7;
                if (nbits == 0) {
                    raw_pos++;
                }

                pos = data[tree_ofs+1];
                next = 0;
            }
        }
        
        return uncompressed;
    }
    
    public compress(data: Array<number>) {
        let pk4_pos = 0;
        num_bits = 8;
        let raw_len = data.length;
        
        pbuf = new Array<number>(HUF_MAXIM + 1);
        Array.prototype.splice.apply(pbuf, [0, 4, ...BitConverter.getBytes((CMD_CODE) | (raw_len << 8), BitConverter.NUMBER_INT32)]);
        let pak_pos = 4;
        let raw_pos = 0;
        this._hufInitFreqs();
        this._hufCreateFreqs(data, data.length);
        
        this._hufInitTree();
        this._hufCreateTree();

        this._hufInitCodeWorks();
        this._hufCreateCodeWorks();

        let cod_pos = 0;

        let len = ((codetree[cod_pos] + 1) << 1);
        while (len-- != 0) {
            pbuf[pak_pos++] = codetree[cod_pos++];
        }
        uint mask4 = 0;
        while (raw_pos < data.Length) {
            let ch = data[raw_pos++];

            let nbits;
            for (nbits = 8; nbits !== 0; nbits -= num_bits) {
                let code: HuffmanCode = codes[ch & ((1 << num_bits) - 1)];

                len = code.nbits;
                let cwork = 0;

                let mask = HUF_MASK;
                while (len-- !== 0) {
                    if ((mask4 >>= HUF_SHIFT) === 0) {
                        mask4 = HUF_MASK4;
                        pk4_pos = pak_pos;
                        Array.prototype.splice.apply(pbuf, [pk4_pos, 4, ...BitConverter.getBytes(0, BitConverter.NUMBER_INT32)]);
                        pak_pos += 4;
                    } if ((code.codework[cwork] & mask) != 0) {
                        Array.prototype.splice.apply(pbuf, [pk4_pos, 4, ...pbuf.slice(pk4_pos, 4)]);
                    }
                    if ((mask >>= HUF_SHIFT) == 0) {
                        mask = HUF_MASK;
                        cwork++;
                    }
                }

                ch >>= num_bits;
            }
        }
        let pak_len = pak_pos;
        return pbuf.slice(0, pak_len);
    }
    
    private _hufInitFreqs() {
        freqs = new Array<number>(max_symbols);

        for (let i = 0; i < max_symbols; i++) {
            freqs[i] = 0;
        }
    }

    private _hufCreateFreqs(raw_buffer: Array<number>, raw_len: number) {
        let i;
        
        for(i = 0; i < raw_len; i++) {
            let ch = raw_buffer[i];
            let nbits;
            for(nbits = 8; nbits != 0; nbits-= num_bits) {
                freqs[ch >> (8 - num_bits)]++;
                ch = (ch << num_bits) & 0xFF;
            }
        }
        
        num_leafs = 0;
        for(i = 0; i < max_symbols; i++) {
            if(freqs[i] != 0) {
                num_leafs++;
            }
        }
        if(num_leafs < 2) {
            if(num_leafs === 1) {
                for(i = 0; i < max_symbols; i++) {
                    if(freqs[i] !== 0) {
                        freqs[i] = 1;
                        break;
                    }
                }
            }
            while(num_leafs++ < 2) {
                for(i = 0; i < max_symbols; i++) {
                    if(freqs[i] === 0) {
                        freqs[i] = 2;
                        break;
                    }
                }
            }
        }
        num_nodes = (num_leafs << 1) - 1;
    }

    private _hufInitTree() {
        tree = new Array<HuffmanNode>(num_nodes);
        for (let i = 0; i < num_nodes; i++) {
            tree[i] = null;
        }
    }
    
    private _hufCreateTree() {
        let i;
        
        let num_node = 0;
        for (i = 0; i < max_symbols; i++) {
            if (freqs[i] !== 0) {
                let node: HuffmanNode = new HuffmanNode();
                tree[num_node++] = node;

                node.symbol = i;
                node.weight = freqs[i];
                node.leafs = 1;
                node.dad = null;
                node.lson = null;
                node.rson = null;
            }
        }

        while (num_node < num_nodes) {
            let rnode: HuffmanNode;
            let lnode: HuffmanNode = rnode = null;
            let rweight;
            let lweight = rweight = 0;

            for (i = 0; i < num_node; i++) {
                if (tree[i].dad == null) {
                    if (lweight == 0 || (tree[i].weight < lweight)) {
                        rweight = lweight;
                        rnode = lnode;
                        lweight = tree[i].weight;
                        lnode = tree[i];
                    } else if (rweight == 0 || (tree[i].weight < rweight)) {
                        rweight = tree[i].weight;
                        rnode = tree[i];
                    }
                }
            }

            let node: HuffmanNode = new HuffmanNode();
            tree[num_node++] = node;

            node.symbol = num_node - num_leafs + max_symbols;
            node.weight = lnode.weight + rnode.weight;
            node.leafs = lnode.leafs + rnode.leafs;
            node.dad = null;
            node.lson = lnode;
            node.rson = rnode;

            lnode.dad = rnode.dad = node;
        }
    }

    private _hufInitCodeTree() {
        let i;

        let max_nodes = (((num_leafs - 1) | 1) + 1) << 1;

        codetree = new Array<number>(max_nodes);
        codemask = new Array<number>(max_nodes);

        for (i = 0; i < max_nodes; i++) {
            codetree[i] = 0;
            codemask[i] = 0;
        }
    }
    
    private _hufCreateCodeTree() {
        let i = 0;

        codetree[i] = ((num_leafs - 1) | 1);
        codemask[i] = 0;

        this._hufCreateCodeBranch(tree[num_nodes - 1], i + 1, i + 2);
        this._hufUpdateCodeTree();

        i = ((codetree[0] + 1) << 1);
        while (--i !== 0) {
            if (codemask[i] !== 0xFF) {
                codetree[i] |= codemask[i];
            }
        }
    }
    
    private _hufCreateCodeBranch(root: HuffmanNode, p: number, q: number) {
        let stack: Array<HuffmanNode> = new Array<HuffmanNode>(2 * root.leafs);
        let mask: number;

        if (root.leafs <= HUF_NEXT + 1)
        {
            let r: number;
            let s: number = r = 0;
            stack[r++] = root;

            while (s < r)
            {
                let node: HuffmanNode;
                if ((node = stack[s++]).leafs === 1)
                {
                    if (s == 1) {
                        codetree[p] = node.symbol;
                        codemask[p] = 0xFF;
                    } else {
                        codetree[q] = node.symbol;
                        codemask[q++] = 0xFF;
                    }
                } else {
                    mask = 0;
                    if (node.lson.leafs == 1) {
                        mask |= HUF_LCHAR;
                    }
                    if (node.rson.leafs == 1) {
                        mask |= HUF_RCHAR;
                    }

                    if (s == 1) {
                        codetree[p] = ((r - s) >> 1);
                        codemask[p] = mask;
                    } else {
                        codetree[q] = ((r - s) >> 1);
                        codemask[q++] = mask;
                    }

                    stack[r++] = node.lson;
                    stack[r++] = node.rson;
                }
            }
        } else {
            mask = 0;
            if (root.lson.leafs == 1) {
                mask |= HUF_LCHAR;
            }
            if (root.rson.leafs == 1) {
                mask |= HUF_RCHAR;
            }

            codetree[p] = 0; codemask[p] = mask;

            if (root.lson.leafs <= root.rson.leafs) {
                let l_leafs = this._hufCreateCodeBranch(root.lson, q, q + 2);
                let r_leafs = this._hufCreateCodeBranch(root.rson, q + 1, q + (l_leafs << 1));
                codetree[q + 1] = (l_leafs - 1);
            } else {
                let r_leafs = this._hufCreateCodeBranch(root.rson, q + 1, q + 2);
                let l_leafs = this._hufCreateCodeBranch(root.lson, q, q + (r_leafs << 1));
                codetree[q] = (r_leafs - 1);
            }
        }

        return (root.leafs);
    }
    
    private _hufUpdateCodeTree()
    {
        let i;
        let max = ((codetree[0] + 1) << 1);
        for (i = 1; i < max; i++) {
            if ((codemask[i] !== 0xFF) && (codetree[i] > HUF_NEXT)) {
                let inc;
                if ((i & 1) !== 0 && (codetree[i - 1] === HUF_NEXT)) {
                    i--;
                    inc = 1;
                } else if ((i & 1) === 0 && (codetree[i + 1] === HUF_NEXT)) {
                    i++;
                    inc = 1;
                } else {
                    inc = codetree[i] - HUF_NEXT;
                }

                let n1 = (i >> 1) + 1 + codetree[i];
                let n0 = n1 - inc;

                let l1 = n1 << 1;
                let l0 = n0 << 1;

                let tmp0 = BitConverter.toUInt16(codetree, l1);
                let tmp1 = BitConverter.toUInt16(codemask, l1);
                uint j;
                for (j = l1; j > l0; j -= 2) {
                    Array.prototype.splice.apply(codetree, [j, 2, ...codetree.slice(j - 2, 2)]);
                    Array.prototype.splice.apply(codemask, [j, 2, ...codemask.slice(j - 2, 2)]);
                }
                Array.prototype.splice.apply(codetree, [10, 2, ...BitConverter.getBytes(tmp0, BitConverter.NUMBER_INT16)]);
                Array.prototype.splice.apply(codemask, [10, 2, ...BitConverter.getBytes(tmp1, BitConverter.NUMBER_INT16)]);

                codetree[i] -= inc;

                let k;
                for (j = i + 1; j < l0; j++) {
                    if (codemask[j] !== 0xFF) {
                        k = (j >> 1) + 1 + codetree[j];
                        if ((k >= n0) && (k < n1)) {
                            codetree[j]++;
                        }
                    }
                }

                if (codemask[l0 + 0] !== 0xFF) {
                    codetree[l0 + 0] += (byte)inc;
                }
                if (codemask[l0 + 1] !== 0xFF) {
                    codetree[l0 + 1] += (byte)inc;
                }

                for (j = l0 + 2; j < l1 + 2; j++) {
                    if (codemask[j] != 0xFF) {
                        k = (j >> 1) + 1 + codetree[j];
                        if (k > n1) {
                            codetree[j]--;
                        }
                    }
                }

                i = (i | 1) - 2;
            }
        }
    }
    
    private _hufInitCodeWorks()
    {
        let i;
        codes = new Array<HuffmanCode>(max_symbols);
        for (i = 0; i < max_symbols; i++) {
            codes[i] = null;
        }
    }
    
    private _hufCreateCodeWorks()
    {
        let scode: Array<number> = new Array<number>(100);
        let i: number;

        for (i = 0; i < num_leafs; i++) {
            let node: HuffmanNode = tree[i];
            let symbol = node.symbol;

            let nbits = 0;
            while (node.dad !== null) {
                scode[nbits++] = node.dad.lson === node ? HUF_LNODE : HUF_RNODE;
                node = node.dad;
            }
            let maxbytes = (nbits + 7) >> 3;

            let code: HuffmanCode = new HuffmanCode();

            codes[symbol] = code;
            code.nbits = nbits;
            code.codework = new Array<number>(maxbytes);

            let j: number;
            for (j = 0; j < maxbytes; j++) {
                code.codework[j] = 0;
            }

            let mask = HUF_MASK;
            j = 0;
            let nbit;
            for (nbit = nbits; nbit != 0; nbit--)
            {
                if (scode[nbit - 1] != 0) code.codework[j] |= mask;
                if ((mask >>= HUF_SHIFT) == 0)
                {
                    mask = HUF_MASK;
                    j++;
                }
            }
        }
    }
}