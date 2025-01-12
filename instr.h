#ifndef INSTR_H
#define INSTR_H

#include <cmath>
#include <cstdint>
#include <map>
#include <string>
class CPU;

class Instr {
private:
  uint32_t m_opcode;
  uint32_t m_r1;
  uint32_t m_r2;
  int32_t m_r3_imm; // sign-extended

public:
  static std::map<std::string, uint32_t> m_instrsDict;
  static void prepareDict() {
#define _ISA(_opcode, _name, _execute, _asmargs, _disasmargs, _dumpregs)       \
  m_instrsDict[#_name] = _opcode;
#include "ISA.h"
#undef _ISA
  }
  void clear();
  void executeCode(CPU *cpu, uint32_t code);
  void execute(CPU *cpu);
  std::string assembler(std::stringstream &input);
  std::string disasm();
  std::string dumpRegs(CPU *cpu);
  uint32_t code();
  void decode(uint32_t code);
  void setOpcode(uint32_t opcode) { m_opcode = opcode; }
};

/* #define _ISA(_opcode, _name, _execute, _asmargs, _disasmargs, _dumpregs)\
class Instr_##_name : public Instr {\
    void execute(CPU& cpu) {_execute};\
    void disasm(std::stringstream &input) {std::string arg; _disasm }\
};
#include "ISA.h"
#undef _ISA */

#endif // INSTR_H
