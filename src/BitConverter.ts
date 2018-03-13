/**
 * Reimplements part of C# BitConverter
 */

export default class BitConverter {
    public static const LITTLE_ENDIAN: number = 0;
    public static const BIG_ENDIAN: number = 1;
    
    public static const NUMBER_AUTO: number = 0;
    public static const NUMBER_BYTE: number = 1;
    public static const NUMBER_INT16: number = 2;
    public static const NUMBER_INT32: number = 3;
    public static const NUMBER_INT64: number = 4;
    
    public static toUint(data: Array<number>, startIndex: number = 0, intLength: number = BitConverter.NUMBER_INT32, endianness: number = BitConverter.LITTLE_ENDIAN) {
        let numberValue = 0;
        let bytesLength = Math.pow(2, actualByteLength - 1);
        for(let idx = 0; idx < bytesLength; idx++) {
            const currentIndex = endianness === BitConverter.BIG_ENDIAN ? startIndex + bytesLength - 1 - idx : startIndex + idx;
            numberValue|= ((data[currentIndex] & 0xFF) << idx);
        }
        return numberValue;
    }
    
    public static toUInt16(data: Array<number>, startIndex: number = 0, endianness: number = BitConverter.LITTLE_ENDIAN) {
        return BitConverter.toUInt(data, startIndex, BitConerter.NUMBER_INT16, endianness);
    }
    
    public static toUInt16LE(data: Array<number>, startIndex: number = 0) {
        return BitConverter.toUInt16(data, startIndex, BitConverter.LITTLE_ENDIAN);
    }
    
    public static toUInt16BE(data: Array<number>, startIndex: number = 0) {
        return BitConverter.toUInt16(data, startIndex, BitConverter.BIG_ENDIAN);
    }
    
    public static toUInt32(data: Array<number>, startIndex: number = 0, endianness: number = BitConverter.LITTLE_ENDIAN) {
        return BitConverter.toUInt(data, startIndex, BitConerter.NUMBER_INT32, endianness);
    }
    
    public static toUInt32LE(data: Array<number>, startIndex: number = 0) {
        return BitConverter.toUInt32(data, startIndex, BitConverter.LITTLE_ENDIAN);
    }
    
    public static toUInt32BE(data: Array<number>, startIndex: number = 0) {
        return BitConverter.toUInt32(data, startIndex, BitConverter.BIG_ENDIAN);
    }
    
    public static getBytes(value: number, byteLength: number = BitConverter.NUMBER_INT32, endianness: number = BitConverter.LITTLE_ENDIAN) {
        let actualByteLength = byteLength;
        if(actualByteLength === 0) {
            let bigVal = 0x100;
            let actualByteLength = 1;
            while(bigVal <= number && actualByteLength < BitConverter.NUMBER_INT64) {
                bigVal*= bigVal;
                actualByteLength++;
            }
        }
        let bytesLength = Math.pow(2, actualByteLength - 1);
        let bytes = new Array<number>(bytesLength);
        
        for(let idx = 0; idx < bytesLength; idx++) {
            const currentIndex = endianness === BitConverter.BIG_ENDIAN ? bytesLength - 1 - idx : idx;
            bytes[currentIndex] = (value >> (idx * 8)) & 0xFF;
        }
        return bytes;
    }
}