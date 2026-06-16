# https://users.ece.cmu.edu/~koopman/crc/#

function koopman_CRC_hex(n::Int, k::Int)

    if k ≤ 0
        throw(error(lazy"k must be greater than zero"))
    end

    if n ≤ k
        throw(error(lazy"n must be greater than k"))
    end

    CRC_size = n - k

    if CRC_size > 32
        throw(error(lazy"CRC size must be less than or equal to 32"))
    elseif CRC_size > 28
        offset = 28

        crc_matrix = [
        #CRC            29         30         31         32 
                0x16dfbf51 0x31342a2f 0x737e312b 0xad0424f3;  
                0x16dfbf51 0x31342a2f 0x737e312b 0xad0424f3;
                0x11c4dfb5 0x2254329d 0x52aa4332 0xc9d204f5;
                0x1cf492f3 0x2adf3aaf 0x74f9e7cb 0xd419cc15;
                0x1cf492f3 0x2ad4a56a 0x74f9e7cb 0x9960034c;
                0x12e8b5b6 0x2a9b3e15 0x60f2920b 0xf8c9140a;
                0x13a46755 0x2017ed6a 0x60f2920b 0xf8c9140a;
                0x1e150a87 0x242c0684 0x6c740b8d 0x9d7f97d6;
                0x1e150a87 0x242c0684 0x6c740b8d 0xb49c1c96;
                0x1c27bd8b 0x34c8e00d 0x456a3501 0x85b9561d;
                0x1c27bd8b 0x2468c69c 0x6bee283f 0x950ebfae;
                0x13a6f65c 0x2b967ef9 0x6624b2eb 0x93b39b1b;
                0x13a6f65c 0x3c9a0b27 0x47e62564 0xa094afb5;
                0x12ff393a 0x290d6d0e 0x52d246e1 0xa2572962;
                0x15e165a6 0x23136e56 0x6d094c5d 0xe89061db;
                0x0        0x229df1ac 0x47d2d9ab 0xa86be4db;
                0x0        0x0        0x46e56a7c 0x973afb51
        ]
        length_matrix = [
        #CRC           29         30         31         32 
      
                536870882 1073741793 2147483616 4294967263;
                268435426  536870881 1073741792 2147483615;
                    16356      32737      32738      65505;
                    16356      16356      32738      32738;
                      484        993        992        992;
                      483        483        992        992;
                      100        100        100        223;
                      100        100        100        100;
                       35         36         36         38;
                       35         35         35         36;
                       14         16         18         20;
                       14         14         16         19;
                        9         11         12         15;
                        8          9         11         13;
                        0          4          5          7;
                        0          0          4          5;
                        0          0          0          0
        ]
    elseif CRC_size > 24
        offset = 24
        crc_matrix = [
        #CRC           25        26        27        28
                0x101690c 0x33c19ef 0x5e04635 0x91dc1e3;
                0x101690c 0x33c19ef 0x5e04635 0x91dc1e3;
                0x10bba2d 0x278b495 0x745e8bf 0xb67b511;
                0x1b9189d 0x2c45446 0x6c3ff0d 0x9037604;
                0x1b9189d 0x2186c30 0x6c3ff0d 0xd120245;
                0x136fd31 0x2bd893b 0x521f64b 0xb9ccb75;
                0x136fd31 0x2bd893b 0x4cb658f 0xb9ccb75;
                0x12b00d4 0x311e9ad 0x4429686 0xeaa72ab;
                0x12b00d4 0x32def69 0x51aff9a 0xacb6aed;
                0x162054b 0x248d3be 0x474fd47 0xb094a3e;
                0x15ed6a9 0x2bfbd8f 0x4258c0f 0xb094a3e;
                0x12728bf 0x2d7a067 0x6986313 0xe9dadcb;
                0x1291ccf 0x23bb612 0x6a611bf 0xaf74fc7;
                0x0       0x251f66b 0x58695e3 0xcf11b95;
                0x0       0x251f66b 0x65bd513 0xcf11b95;
        ]
        length_matrix = [
        #CRC      25        26        27        28
            33554406  67108837 134217700  268435427;
            16777190  33554405  67108836  134217699;
                4072      8165      8166      16357;
                4072      4072      8166       8166;
                 230       230       484        483;
                 230       230       230        483;
                  40        41        48         99;
                  40        40        41         48;
                  24        24        36         35;
                  23        23        23         35;
                   8         9        11         15;
                   7         8         9         11;
                   0         5         7          8;
                   0         5         6          8;
                   0         0         0          0;
        ] 
    elseif CRC_size > 16   
        offset = 16
        crc_matrix = [
        #CRC         17      18      19      20       21       22       23       24
                0x16fa7 0x23979 0x6fb57 0xb5827 0x1707ea 0x308fd3 0x540df0 0x8f90e3;
                0x16fa7 0x23979 0x6fb57 0xb5827 0x1707ea 0x308fd3 0x540df0 0x8f90e3;
                0x1165d 0x25f53 0x77b0f 0xc1acf 0x10df8f 0x248794 0x400154 0x9945b1;
                0x1724e 0x39553 0x5685a 0xc8a89 0x1edfb7 0x2a952a 0x6bc0f5 0x98ff8c;
                0x1724e 0x32c69 0x5685a 0xe2023 0x1edfb7 0x395b53 0x6bc0f5 0xbd80de;
                0x1751b 0x25f6a 0x50b49 0x8810e 0x12faa5 0x289cfe 0x5e2419 0x880ee6;
                0x11bf5 0x25f6a 0x779c7 0xd41cf 0x198313 0x289cfe 0x469d7c 0xcba785;
                0x123bd 0x27bbc 0x7573f 0xbe73e 0x16e976 0x2aedd3 0x53df6e 0xed93bb;
                0x176a7 0x2e7de 0x44f75 0xe6233 0x16e976 0x247bc4 0x463b77 0xc7ad89;
                0x0     0x26a3d 0x6d133 0x8d3cc 0x165751 0x36f627 0x463b77 0x8cd929;
                0x0     0x0     0x51d79 0x9d587 0x165751 0x22efb7 0x49ad52 0x8cd929;
                0x0     0x0     0x0     0x0     0x0      0x25d467 0x4b79d1 0xd9588b;
                0x0     0x0     0x0     0x0     0x0      0x0      0x4b79d1 0xb73e91;
        ]
        length_matrix = [
        #CRC    17     18     19      20      21      22      23       24
            131054 262125 524268 1048555 2097130 4194281 8388584 16777191;
             65518 131053 262124  524267 1048554 2097129 4194280  8388583;
               240    493    494    1005    1004    2025    2026     4073;
               240    240    494     494    1004    1004    2026     2026;
                46     45     46      49     106     105     106      231;
                22     45     45      45      48     105     105      105;
                 8     11     13      21      20      22      26       39;
                 6      8     10      13      20      20      24       26;
                 0      5      7      11      10      12      24       23;
                 0      0      5       7      10      10      13       23;
                 0      0      0       0       0       5       5        7;
                 0      0      0       0       0       0       5        6;
                 0      0      0       0       0       0       0        0;
        ]   
    elseif CRC_size > 2
        offset = 2
        crc_matrix = [
        #CRC      3    4     5    6    7    8      9    10    11    12     13     14     15     16
                0x5  0x9  0x12 0x33 0x65 0xe7  0x119 0x327 0x5db 0x987 0x1abf 0x27cf 0x4f23 0x8d95;
                0x5  0x9  0x12 0x33 0x65 0xe7  0x119 0x327 0x5db 0x987 0x1abf 0x27cf 0x4f23 0x8d95;
                0x0  0x0  0x15 0x23 0x5b 0x98  0x17d 0x247 0x583 0x8f3 0x12e6 0x2322 0x4306 0xd175;
                0x0  0x0  0x0  0x0  0x72 0xeb  0x185 0x2b9 0x5d7 0xbae 0x1e97 0x212d 0x6a8d 0xac9a;
                0x0  0x0  0x0  0x0  0x0  0x9b  0x13c 0x28e 0x532 0xb41 0x1e97 0x372b 0x573a 0x9eb2;
                0x0  0x0  0x0  0x0  0x0  0x0   0x0   0x29b 0x571 0xa4f 0x12a5 0x28a9 0x5bd5 0x968b;
                0x0  0x0  0x0  0x0  0x0  0x0   0x0   0x0   0x4f5 0xa4f 0x10b7 0x2371 0x630b 0x8fdb;
                0x0  0x0  0x0  0x0  0x0  0x0   0x0   0x0   0x0   0x0   0x0    0x0    0x5a47 0xe92f;
                0x0  0x0  0x0  0x0  0x0  0x0   0x0   0x0   0x0   0x0   0x0    0x0    0x0    0xed2f;
        ]
        length_matrix = [
        #CRC 3     4     5     6      7     8     9     10    11    12    13     14    15    16
            4    11    26    57    120   247   502   1013  2036  4083  8178  16369 32752 65519;
            0     0    10    25     56   119   246    501  1012  2035  4082   8177 16368 32751;
            0     0     0     0      4     9    13     21    26    53    52    113   136   241;
            0     0     0     0      0     4     8     12    22    27    52     57   114   135;
            0     0     0     0      0     0     0      5    12    11    12     13    16    19;
            0     0     0     0      0     0     0      0     4    11    11     11    12    15;
            0     0     0     0      0     0     0      0     0     0     0      0     5     6;
            0     0     0     0      0     0     0      0     0     0     0      0     0     5;
            0     0     0     0      0     0     0      0     0     0     0      0     0     0;
        ]    
    else
        throw(error(lazy"CRC size must be greater than 2"))
    end

    hd, i = find_hd(k,length_matrix,CRC_size,offset)

    return hd, crc_matrix[i,CRC_size-offset]

end

function find_hd(
    k::Int,
    length_matrix::Matrix{Int},
    CRC_size::Int,
    offset::Int)

    i = 1
    j = CRC_size-offset
    max_length = length_matrix[i,j]

    if k > max_length
        return 2, 1
    else
        while max_length > 0
            i +=1
            max_length = length_matrix[i,j]
            if k > max_length
                return i + 1, i
            end
        end
        return i + 1, i
    end
end