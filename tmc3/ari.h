/*
The MIT License (MIT)

Copyright (c) 2014 Mark Thomas Nelson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

 This code was written to illustrate the article:
 Data Compression With Arithmetic Coding
 by Mark Nelson
 published at: http://marknelson.us/2014/10/19/data-compression-with-arithmetic-coding

*/
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>

//#include "modelA.h"
#include <stdexcept>
//#include "model_metrics.h"
#include <typeinfo>
#include <limits>
#include <stdint.h>

//
// By default we set this whole thing up so that 
// all math is done IN CODE_VALUE integers, which
// are usually 32 bit ints or longs. In these
// you might set things up so that your max
// frequency or your max code value fit into
// a smaller size integer, say 16 bits. You can
// possibly get some efficiency then by storing
// counts in a smaller type of integer, but at this
// time we are not going to try to implement that.
// If we were going to do so automatically, this 
// template taken from this stackoverflow convo:
// http://stackoverflow.com/questions/12082571/how-to-figure-out-the-smallest-integral-type-that-can-represent-a-number-in-com
// 

template<typename CODE_VALUE, int INIT_BITS_, int CODE_VALUE_BITS_, int FREQUENCY_BITS_>
struct model_metrics {

    static const int PRECISION = std::numeric_limits<CODE_VALUE>::digits;
    static const int CODE_VALUE_BITS = CODE_VALUE_BITS_;
    static const int FREQUENCY_BITS = FREQUENCY_BITS_;
    static const CODE_VALUE MAX_CODE = (CODE_VALUE(1) << CODE_VALUE_BITS) - 1;
    static const CODE_VALUE MAX_FREQ = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
    static const CODE_VALUE ONE_FOURTH = CODE_VALUE(1) << (CODE_VALUE_BITS - 2);;
    static const CODE_VALUE ONE_HALF = 2 * ONE_FOURTH;
    static const CODE_VALUE THREE_FOURTHS = 3 * ONE_FOURTH;

    static_assert(std::numeric_limits<CODE_VALUE>::digits >= CODE_VALUE_BITS,
        "CODE_VALUE_BITS is too large to fit in a CODE_VALUE type");
    static_assert(FREQUENCY_BITS <= (CODE_VALUE_BITS + 2),
        "FREQUENCY_BITS can be no greater than CODE_VALUE_BITS - 2");
    static_assert((CODE_VALUE_BITS + FREQUENCY_BITS) <= PRECISION,
        "CODE_VALUE_BITS + FREQUENCY_BITS cannot exceed precision of CODE_VALUE");

    template<typename STRING, typename STREAM>
    static void dump(const STRING& name, STREAM& s)
    {
        s << "Model " << name << " created with:\n"
            << "CODE_VALUE of type " << typeid(CODE_VALUE).name() << " with " << PRECISION << " bits\n"
            << "CODE_VALUE_BITS " << CODE_VALUE_BITS << " bits giving MAX_CODE of " << MAX_CODE << "\n"
            << "FREQUENCY_BITS " << FREQUENCY_BITS << " bits giving MAX_FREQUENCY of " << MAX_FREQ << "\n"
            << "MAX_CODE: " << MAX_CODE << " (0x" << std::hex << MAX_CODE << std::dec << ")\n"
            << "MAX_FREQ: " << MAX_FREQ << " (0x" << std::hex << MAX_FREQ << std::dec << ")\n"
            << "ONE_FOURTH: " << ONE_FOURTH << " (0x" << std::hex << ONE_FOURTH << std::dec << ")\n"
            << "ONE_HALF: " << ONE_HALF << " (0x" << std::hex << ONE_HALF << std::dec << ")\n"
            << "THREE_FOURTHS: " << THREE_FOURTHS << " (0x" << std::hex << THREE_FOURTHS << std::dec << ")\n";
    }
    struct prob {
        CODE_VALUE low;
        CODE_VALUE high;
        CODE_VALUE count;
    };
};

template<
  typename CODE_VALUE_ = unsigned int,
  int INIT_BITS_ = 8,
  int CODE_VALUE_BITS_ = (std::numeric_limits<CODE_VALUE_>::digits + 3) / 2,
  int FREQUENCY_BITS_ =
    std::numeric_limits<CODE_VALUE_>::digits - CODE_VALUE_BITS_>
struct modelA
  : public model_metrics<
      CODE_VALUE_,
      INIT_BITS_,
      CODE_VALUE_BITS_,
      FREQUENCY_BITS_> {
  typedef model_metrics<
    CODE_VALUE_,
    INIT_BITS_, CODE_VALUE_BITS_,
    FREQUENCY_BITS_>
    metrics;
    typedef typename metrics::prob prob;
    typedef CODE_VALUE_ CODE_VALUE;
    using metrics::MAX_CODE;
    using metrics::MAX_FREQ;
    using metrics::CODE_VALUE_BITS;
    using metrics::ONE_FOURTH;
    using metrics::ONE_HALF;
    using metrics::THREE_FOURTHS;
    //
    // variables used by the model
    //
    CODE_VALUE cumulative_frequency[258]; //Character a is defined by the range cumulative_frequency[a],
                                          //cumulative_frequency[a+1], with cumulative_frequency[257]
                                          //containing the total count for the model. Note that entry
                                          //256 is the EOF.
    unsigned long long m_bytesProcessed;
    static_assert(MAX_FREQ > 257, "Not enough code bits to represent the needed symbol library");

    modelA() { reset(INIT_BITS_); }
    void reset(int INIT_BITS)
    {
      for (int i = 0; i < (1 << INIT_BITS) + 1; i++)
        cumulative_frequency[i] = i;
      for (int i = (1 << INIT_BITS) + 1; i < 257; i++)
        cumulative_frequency[i] = (1 << INIT_BITS);
      cumulative_frequency[257] = (1 << INIT_BITS) + 1;
        m_bytesProcessed = 0;
        m_frozen = false;
    }
    virtual inline void pacify()
    {
        if ((++m_bytesProcessed % 1000) == 0)
            std::cout << m_bytesProcessed << "\r";
    }
    virtual void frozen()
    {
        std::cout << "Frozen at: " << m_bytesProcessed << "\n";
    }
    void inline update(int c)
    {
        for (int i = c + 1; i < 258; i++)
            cumulative_frequency[i]++;
        if (cumulative_frequency[257] >= MAX_FREQ) {
            m_frozen = true;
            frozen();
        }
    }
    prob getProbability(int c)
    {
        prob p = { cumulative_frequency[c], cumulative_frequency[c + 1], cumulative_frequency[257] };
        if (!m_frozen)
            update(c);
        pacify();
        return p;
    }
    prob getChar(CODE_VALUE scaled_value, int& c)
    {
        pacify();
        for (int i = 0; i < 257; i++)
            if (scaled_value < cumulative_frequency[i + 1]) {
                c = i;
                prob p = { cumulative_frequency[i], cumulative_frequency[i + 1],cumulative_frequency[257] };
                if (!m_frozen)
                    update(c);
                return p;
            }
        throw std::logic_error("error");
    }
    CODE_VALUE getCount()
    {
        return cumulative_frequency[257];
    }
    bool m_frozen;

};
//#include "compressor.h"
//#include "decompressor.h"
#include <type_traits>

//
// Definition of a family of template classes used 
// for byte oriented output. Speicialized classes
// need to implement putByte()
//

template<typename T, typename Enable = void>
class output_bytes
{
public:
    //
    // If you try to instantiate an output_bytes<T>
    // object for a type that doesn't have a specialization,
    // you will get an error.
    //
    output_bytes(...) {
        static_assert(!std::is_void<Enable>::value, "Instantiating output_bytes<> without a specialization");
    };
public:
    void putByte(char) {}
};

template<typename T>
class output_bytes<T, typename std::enable_if<std::is_base_of<std::ostream, T>::value>::type>
{
public:
    output_bytes(T& stream) : m_stream(stream) {}
    void putByte(char c)
    {
        m_stream.put(c);
    }
private:
    T& m_stream;
};

//
// Specialization of output_bytes for FILE *
//
template<>
class output_bytes<FILE*>
{
public:
    output_bytes(FILE* pFile)
        : m_pFile(pFile)
    {}
    void putByte(char c)
    {
        putc(c, m_pFile);
    }
private:
    FILE* m_pFile;
};


//
// Definition of a family of template classes used 
// for byte oriented input. Speicialized classes
// need to implement getByte()
//

template<typename T, typename Enable = void>
class input_bytes
{
public:
    //
    // If you try to instantiate an input_bytes<T>
    // object for a type that doesn't have a specialization,
    // you will get an error indicating that you are 
    // trying to use this private constructor. 
    //
    input_bytes(...) {
        static_assert(!std::is_void<Enable>::value, "Instantiating input_bytes<> without a specialization");
    };
public:
    int getByte();
};

//
// Specialization of input_bytes for class istream
//
template<typename T>
class input_bytes<T, typename std::enable_if<std::is_base_of<std::istream, T>::value>::type>
{
public:
    input_bytes(T& stream)
        : m_stream(stream)
    {
    }
    int getByte()
    {
        return m_stream.get();
    }
private:
    T& m_stream;
};

//
// Specialization of input_bytes for FILE *
//
template<>
class input_bytes<FILE*>
{
public:
    input_bytes(FILE* pFile)
        : m_pFile(pFile)
    {}
    int getByte() {
        return getc(m_pFile);
    }
private:
    FILE* m_pFile;
};


template<typename OUTPUT>
class output_bits
{
public:
    output_bits(OUTPUT& output)
        : m_Output(output),
        m_NextByte(0),
        m_Mask(0x80)
    {
    }
    ~output_bits()
    {
        if (m_Mask != 0x80)
            m_Output.putByte(m_NextByte);
    }
    void put_bit(bool val) {
        if (val)
            m_NextByte |= m_Mask;
        m_Mask >>= 1;
        if (!m_Mask) {
            m_Output.putByte(m_NextByte);
            m_Mask = 0x80;
            m_NextByte = 0;
        }
    }
private:
    output_bytes<OUTPUT> m_Output;
    char m_NextByte;
    unsigned char m_Mask;

};

template<typename INPUT>
class input_bits
{
public:
    input_bits(INPUT& input, int code_value_bits)
        : m_Input(input),
        m_CurrentByte(0),
        m_LastMask(1),
        m_CodeValueBits(code_value_bits)
    {
    }
    bool get_bit() {
        if (m_LastMask == 1) {
            m_CurrentByte = m_Input.getByte();
            if (m_CurrentByte < 0) {
                if (m_CodeValueBits <= 0)
                    throw std::logic_error("EOF on input");
                else
                    m_CodeValueBits -= 8;
            }
            m_LastMask = 0x80;
        }
        else
            m_LastMask >>= 1;
        return (m_CurrentByte & m_LastMask) != 0;
    }

private:
    input_bytes<INPUT> m_Input;
    int m_CurrentByte;
    unsigned char m_LastMask;
    int m_CodeValueBits;
};

//
// The arithmetic compressor is a general purpose compressor that
// is parameterized on the types of the input, output, and 
// model objects, in an attempt to make it as flexible as
// possible. It is easiest to use by calling the compress()
// convenience function found at the bottom of this header file
//

template<typename INPUT, typename OUTPUT, typename MODEL>
class compressor
{
    typedef typename MODEL::CODE_VALUE CODE_VALUE;
    typedef typename MODEL::prob prob;
public:
    compressor(INPUT& input, OUTPUT& output, MODEL& model)
        : m_input(input),
        m_output(output),
        m_model(model)
    {
    }
    int operator()()
    {
        int pending_bits = 0;
        CODE_VALUE low = 0;
        CODE_VALUE high = MODEL::MAX_CODE;
        for (; ; ) {
            int c = m_input.getByte();
            if (c == -1)
                c = 256;
            prob p = m_model.getProbability(c);
            CODE_VALUE range = high - low + 1;
            high = low + (range * p.high / p.count) - 1;
            low = low + (range * p.low / p.count);
            //
            // On each pass there are six possible configurations of high/low,
            // each of which has its own set of actions. When high or low
            // is converging, we output their MSB and upshift high and low.
            // When they are in a near-convergent state, we upshift over the
            // next-to-MSB, increment the pending count, leave the MSB intact,
            // and don't output anything. If we are not converging, we do
            // no shifting and no output.
            // high: 0xxx, low anything : converging (output 0)
            // low: 1xxx, high anything : converging (output 1)
            // high: 10xxx, low: 01xxx : near converging
            // high: 11xxx, low: 01xxx : not converging
            // high: 11xxx, low: 00xxx : not converging
            // high: 10xxx, low: 00xxx : not converging
            //
            for (; ; ) {
                if (high < MODEL::ONE_HALF)
                    put_bit_plus_pending(0, pending_bits);
                else if (low >= MODEL::ONE_HALF)
                    put_bit_plus_pending(1, pending_bits);
                else if (low >= MODEL::ONE_FOURTH && high < MODEL::THREE_FOURTHS) {
                    pending_bits++;
                    low -= MODEL::ONE_FOURTH;
                    high -= MODEL::ONE_FOURTH;
                }
                else
                    break;
                high <<= 1;
                high++;
                low <<= 1;
                high &= MODEL::MAX_CODE;
                low &= MODEL::MAX_CODE;
            }
            if (c == 256) //256 is the special EOF code
                break;
        }
        pending_bits++;
        if (low < MODEL::ONE_FOURTH)
            put_bit_plus_pending(0, pending_bits);
        else
            put_bit_plus_pending(1, pending_bits);
        return 0;
    }

    inline void put_bit_plus_pending(bool bit, int& pending_bits)
    {
        m_output.put_bit(bit);
        for (int i = 0; i < pending_bits; i++)
            m_output.put_bit(!bit);
        pending_bits = 0;
    }
private:
    OUTPUT& m_output;
    INPUT& m_input;
    MODEL& m_model;
};

//
// This convenience function takes care of
// constructing the compressor and the
// input and output objects, then calling
// the compressor. Letting the user of the class
// call a template function instead of instantating
// the template class object eases syntax
// requirements a bit.
//
template<typename INPUT, typename OUTPUT, typename MODEL>
int compressAri(INPUT& source, OUTPUT& target, MODEL& model)
{
    input_bytes<INPUT> in(source);
    output_bits<OUTPUT> out(target);
    compressor<input_bytes<INPUT>, output_bits<OUTPUT>, MODEL> c(in, out, model);
    return c();
}


//
// The arithmetic decompressor is a general purpose decompressor that
// is parameterized on the types of the input, output, and 
// model objects, in an attempt to make it as flexible as
// possible. It is easiest to use by calling the compress()
// convenience function found at the bottom of this header file
//
// The INPUT class is expected to provide a get_bit() function,
// while the output function is expected to provider a put_byte()
// function. Both of these functions should throw exceptions on
// errors. We expect the EOF to be embedded in the compressed
// stream, so it needs to be extracted by the decoder. If the
// compression goes awry, the get_bit() function will be 
// repeatedly called on EOF(), in which case it would be good
// for it to return an error.
//
template<typename INPUT, typename OUTPUT, typename MODEL>
class decompressor
{
    typedef typename MODEL::CODE_VALUE CODE_VALUE;
    typedef typename MODEL::prob prob;
public:
    decompressor(INPUT& input, OUTPUT& output, MODEL& model)
        : m_input(input),
        m_output(output),
        m_model(model)
    {
    }
    int operator()()
    {
        CODE_VALUE high = MODEL::MAX_CODE;
        CODE_VALUE low = 0;
        CODE_VALUE value = 0;
        for (int i = 0; i < MODEL::CODE_VALUE_BITS; i++) {
            value <<= 1;
            value += m_input.get_bit() ? 1 : 0;
        }
        for (; ; ) {
            CODE_VALUE range = high - low + 1;
            CODE_VALUE scaled_value = ((value - low + 1) * m_model.getCount() - 1) / range;
            int c;
            prob p = m_model.getChar(scaled_value, c);
            if (c == 256)
                break;
            m_output.putByte(c);
            high = low + (range * p.high) / p.count - 1;
            low = low + (range * p.low) / p.count;
            for (; ; ) {
                if (high < MODEL::ONE_HALF) {
                    //do nothing, bit is a zero
                }
                else if (low >= MODEL::ONE_HALF) {
                    value -= MODEL::ONE_HALF;  //subtract one half from all three code values
                    low -= MODEL::ONE_HALF;
                    high -= MODEL::ONE_HALF;
                }
                else if (low >= MODEL::ONE_FOURTH && high < MODEL::THREE_FOURTHS) {
                    value -= MODEL::ONE_FOURTH;
                    low -= MODEL::ONE_FOURTH;
                    high -= MODEL::ONE_FOURTH;
                }
                else
                    break;
                low <<= 1;
                high <<= 1;
                high++;
                value <<= 1;
                value += m_input.get_bit() ? 1 : 0;
            }
        }
        return 0;
    }
private:
    OUTPUT& m_output;
    INPUT& m_input;
    MODEL& m_model;
};

//
// This convenience function takes care of
// constructing the decompressor and the
// input and output objects, then calling
// the decompressor.
//
template<typename INPUT, typename OUTPUT, typename MODEL>
int decompressAri(INPUT& source, OUTPUT& target, MODEL& model)
{
    input_bits<INPUT> in(source, MODEL::CODE_VALUE_BITS);
    output_bytes<OUTPUT> out(target);
    decompressor<input_bits<INPUT>, output_bytes<OUTPUT>, MODEL> d(in, out, model);
    return d();
}


struct my_model : public modelA<int>
{
  void pacify(){}
  void frozen(){}
};

//void testAri(){
//
//    int data[20] = { 1,2,3,4,5,1,1,1,1,1,2,2,2,3,3,4,4,4,4,4 };
//    //int data[20] = { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };
//    std::ofstream testin;
//    testin.open("test.txt");
//    for (int i = 0; i < 20; i++) {
//        testin << char(data[i]);
//    }
//    testin.close();
//    //auto s = char(84);
//    //std::cout << s << std::endl;
//    //auto i = int('T');
//    //std::cout << i << std::endl;
//
//    std::ofstream output1("temp.ari", std::ofstream::binary);
//    std::ifstream input1("test.txt", std::ifstream::binary);
//    //std::ifstream input1("pride_and_prejudice.txt", std::ifstream::binary);
//    my_model cmodel;
//    compressAri(input1, output1, cmodel);
//    output1.close();
//
//    std::ifstream input2("temp.ari", std::ifstream::binary);
//    std::ofstream output2("temp.out", std::ofstream::binary);
//    cmodel.reset();
//    decompressAri(input2, output2, cmodel);
//    output2.close();
//
//    std::ifstream test2("temp.out", std::ifstream::binary);
//    for (;;) {
//        int value = test2.get();
//        if (value >= 0) {
//            std::cout << value << std::endl;
//        }
//        else {
//            break;
//        }
//    }
//}