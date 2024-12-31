from collections import Counter, defaultdict
import re
from loguru import logger
from .lipid_class import lipid_class_dict
from .common_molecular import *
from .ion import *

class FormulaParser:

    def __init__(self):
        self.atom_list = None
        self.mass = 0
        self.formula = None

    def split_lipid_string(self, s):
        """
        将最原始的脂质名称分割为多个部分，例如TG (O-50:1) [NL-16:0]分割为['TG (O-50:1)','[NL-16:0]']
        :param s: 最原始的脂质名称，例如TG (O-50:1) [NL-16:0]
        :return:
        """
        parts = []
        stack = []
        last_index = 0

        for i, char in enumerate(s):
            if char in '([{':
                stack.append((char, i))
            elif char in ')]}':
                if stack:
                    open_char, open_index = stack.pop()
                    if open_char == '(' and char == ')':
                        if not stack:
                            parts.append(s[last_index:open_index])
                            parts.append(s[open_index:i + 1])
                            last_index = i + 1
                    elif open_char == '[' and char == ']':
                        if not stack:
                            parts.append(s[last_index:open_index])
                            parts.append(s[open_index:i + 1])
                            last_index = i + 1

        if last_index < len(s):
            parts.append(s[last_index:])

        parts = [part.strip() for part in parts if part.strip()]

        i = 0
        while i < len(parts) - 1:
            if '(' not in parts[i]:
                if '(' not in parts[i + 1]:
                    parts[i + 1] = '(' + parts[i + 1] + ')'
                parts[i] += parts[i + 1]
                del parts[i + 1]
            else:
                i += 1

        return parts

    def merge_element_counts(self, *args):
        result = defaultdict(int)

        for arg in args:
            if isinstance(arg, dict):
                for key, value in arg.items():
                    result[key] += value
            elif isinstance(arg, list):
                for dictionary in arg:
                    for key, value in dictionary.items():
                        result[key] += value
            else:
                raise ValueError(
                    "Unsupported argument type. Each argument must be a dictionary or a list of dictionaries.")

        return dict(result)

    def delete_element_counts(self, main_dict, *args):
        """
        从main_dict中减去args中的元素数目
        :param main_dict:主字典，其中包含元素及其数量的键值对
        :param args:可变数量的参数，每个参数可以是一个字典或包含多个字典的列表。每个字典都包含希望从 main_dict 中减去的元素数量
        :return:更新后的 main_dict，已经减去了 args 中指定的元素数目
        """
        for arg in args:
            if isinstance(arg, dict):
                for key, value in arg.items():
                    main_dict[key] = main_dict.get(key, 0) - value
                    if main_dict[key] <= 0:
                        del main_dict[key]
            elif isinstance(arg, list):
                for dictionary in arg:
                    for key, value in dictionary.items():
                        main_dict[key] = main_dict.get(key, 0) - value
                        if main_dict[key] <= 0:
                            del main_dict[key]
            else:
                raise ValueError("数据格式有误")

        return main_dict

    def remove_outer_brackets(self, s):
        """
        去除字符串两端的括号
        :param s: 脂质名称分割后的结果，例如(18:0/0:0/0:0)、(d5)
        :return: 去除括号后的字符串
        """

        if (s.startswith('(') and s.endswith(')')) or (s.startswith('（') and s.endswith('）')):
            return s[1:-1]
        return s

    def detect_lipid_category(self, lipid_string, lipid_class = lipid_class_dict):
        """
        检测脂质的类别，例如TG、DG、MG、PC、PE等
        :param lipid_string: 脂质名称，根据名称首字母匹配脂质类别，例如TG (O-50:1) [NL-16:0]
        :param lipid_class: 来源于LipidMAPS的脂质类别字典，（Category - Main Class）->（MainClass - SubClass）
        :return: 返回具体的脂质类别，返回的内容为(MainClass, SubClass, Formula)
        """

        match = re.match(r"(\w+)", lipid_string)
        if match:
            lipid_abbreviation = match.group(1)
        else:
            # If no abbreviation is detected, return None
            return None, None, lipid_string

        # Check each category in the lipid_class dictionary
        for category, abbreviations in lipid_class.items():
            if lipid_abbreviation in abbreviations:
                return category, lipid_abbreviation, lipid_string

        # Return None if no category is found
        return None, None, lipid_string

    def check_chain(self, chain):

        """
        检查链的格式是否正确，并返回链的长度、饱和度、缩醛标记、甲基化标记

        :param chain: 单个脂肪酸链，例如O-14:0(13Me)、20:5(5Z,8Z,11Z,14Z,17Z)
        :return: 长度:饱和度、缩醛标记（O:1，P：2）、额外甲基数量（Me:1）
        :rtype: list
        """
        logger.debug(f"{self.formula} - Checking chain: {chain}")

        # 检查链的长度和饱和度
        rc1 = re.search(r'(\d+):(\d+)', chain).group()

        # 确定缩醛标记"O-"、"A-"前者返回1，后者返回2
        rc2 = re.search(r'[OaPp]-', chain)
        if rc2 is None:
            rc2 = 0
        else:
            rc2 = rc2.group()
            if rc2 == "O-" or rc2 == "a-":
                rc2 = 1
            elif rc2 == "P-" or rc2 == "p-":
                rc2 = 2

        # 确定甲基化标记"Me"，返回甲基数量
        rc3 = len(re.findall(r'Me', chain))

        # 确定鞘脂骨架的羟基标记
        rc4 = re.fullmatch(r'([mdt])(\d+):(\d+).*', chain)
        if rc4:
            rc4 = rc4.group(1)
        else:
            rc4 = None

        # 确定脂肪酸链的羟基标记，返回羟基数量
        rc5 = len(re.findall(r'OH', chain))

        return [rc1, rc2, rc3, rc4, rc5]

    def calculate_fatty_acid_chain(self, chain_list: str | list, actual_chain_num: int = None) -> list:
        """
        计算脂肪酸链的碳、氢、氧原子数目，如果输入单个链，则返回单个链的结果，如果输入多个链，则返回多个链的结果。
        :param chain_list: 脂肪酸链的简写形式，支持单个或多个链的输入，例如"16:0"或["16:0", "18:1"]
        :param lipid_sub_class: 脂质子类别，例如TG、DG、MG等，用于校正多链合一的情况下总氧原子数目
        :return: 脂肪酸链的碳、氢、氧原子数目
        """

        if isinstance(chain_list, str):
            chain_list = [chain_list]

        output_chain_list = []
        # 保证每条链均有效
        for chain in [chain for chain in chain_list if chain != '0:0']:

            # 分解输入字符串，获取碳和双键的数量
            carbon_atoms, double_bonds = map(int, chain.split(':'))

            # 根据公式计算氢原子的数量
            hydrogen_atoms = 2 * carbon_atoms - 2 * double_bonds

            output_chain_list.append({'C': carbon_atoms, 'H': hydrogen_atoms, 'O':2 })

        # 根据真实的链数目，校正氧原子数目
        if actual_chain_num is not None and len(chain_list) == 1:
            output_chain_list[0]['O'] += 2 * (actual_chain_num - len(chain_list))


        # 如果是TG或DG，需要补充氧原子数目
        # if lipid_sub_class == "TG" and len(chain_list) == 1:
        #     output_chain_list[0]['O'] += 2 * (3 - 1)
        # elif lipid_sub_class == "DG" and len(chain_list) == 1:
        #     output_chain_list[0]['O'] += 2 * (2 - 1)

        return output_chain_list

    def add_lipid_chains(self, chain1, chain2):
        """
        用于计算两个脂肪酸链的碳、双键之和
        :param chain1: 脂肪酸链1，例如"16:0"
        :param chain2: 脂肪酸链2，例如"18:1"
        :return: 两个脂肪酸链的碳、双键之和，例如"34:1"
        """
        # 解析第一个链的长度和不饱和度
        length1, unsaturation1 = map(int, chain1.split(':'))
        # 解析第二个链的长度和不饱和度
        length2, unsaturation2 = map(int, chain2.split(':'))

        # 将相应的长度和不饱和度进行相加
        total_length = length1 + length2
        total_unsaturation = unsaturation1 + unsaturation2

        # 构造结果字符串
        result = f"{total_length}:{total_unsaturation}"

        return result

    def delete_lipid_chains(self, chain1, chain2):
        """
        用于计算两个脂肪酸链的碳、双键之差
        :param chain1: 脂肪酸链1，例如"34:1"
        :param chain2: 脂肪酸链2，例如"18:1"
        :return: 两个脂肪酸链的碳、双键之和，例如"16:0"
        """
        # 解析第一个链的长度和不饱和度
        length1, unsaturation1 = map(int, chain1.split(':'))
        # 解析第二个链的长度和不饱和度
        length2, unsaturation2 = map(int, chain2.split(':'))

        # 将相应的长度和不饱和度进行相加
        total_length = length1 - length2
        total_unsaturation = unsaturation1 - unsaturation2

        if total_length < 0 or total_unsaturation < 0:
            logger.warning(f"{chain1} - {chain2} = {total_length}:{total_unsaturation}")

        # 构造结果字符串
        result = f"{total_length}:{total_unsaturation}"

        return result

    def split_lipid(self,formula_list: tuple | list):

        """
        将整体的脂肪酸链分割为单个的脂肪酸链，以碳原子数:双键数目的形式返回
        :param formula_list: [('GL', 'TG', '(16:1(9Z)/16:1(9Z)/17:1(9Z))')]
        :return:["16:1","16:1","17:1"]
        """

        if isinstance(formula_list, tuple):
            formula_list = [formula_list]

        res = []

        for formula in formula_list:
            bone = formula[1]
            formula_str = formula[2]
            logger.debug(f"{self.formula} - formula_str: {formula_str}")

            if formula_str.startswith(bone):
                chains = self.remove_outer_brackets(formula_str[len(bone):])
                logger.debug(f"{self.formula} - Total chains list: {chains}")
                chain_list = re.split(r'[/_]', chains)
                chain_list = [self.check_chain(chain) for chain in chain_list]
                res.append((formula[0], formula[1], formula[2], chain_list))
            else:
                logger.error(f"{self.formula} - ERROR: {formula_str} is not a valid lipid formula")

        return res

    ##
    def countDeuterium(self, formula_list: str | list):
        """
        计算初次分割之后，脂质中氘原子数目
        :param formula_list: ["TG(20:5/22:6/20:5)","(d5)"]
        :return: {'D':5}
        """
        if isinstance(formula_list, str):
            formula_list = [formula_list]

        deuterium_counts = []

        for formula in formula_list:
            matches = re.finditer(r'\bd(\d+)\b(?![:/])', formula)
            for match in matches:
                number = match.group(1)
                deuterium_counts.append(int(number))  # Convert string to int and add to list

        return {'D': sum(deuterium_counts)} if deuterium_counts else {}

    def countNL(self, formula_list: str | list):
        """
        查找并返回初次分割之后脂质中的中性丢失链
        :param formula_list: ["TG (O-50:1)","[NL-16:0]"]
        :return: ["16:0"]
        """
        # Ensure input is a list
        if isinstance(formula_list, str):
            formula_list = [formula_list]

        nl_chains = []
        for formula in formula_list:
            match = re.search(r'[\[\(]NL-(\d+:\d+)[\]\)]', formula)
            if match:
                nl_chains.append(match.group(1))

        nl_chains = self.calculate_fatty_acid_chain(nl_chains)
        total_nl_chain = self.merge_element_counts(nl_chains)

        return total_nl_chain

    def countIon(self, formula_list: str | list):
        """
        查找并返回初次分割之后脂质中的加和离子
        :param formula_list:
        :return:
        """
        if isinstance(formula_list, str):
            formula_list = [formula_list]

        matched_ions = [ion_list[ion] for ion in formula_list if ion in ion_list.keys()]

        if len(matched_ions) == 1:
            return matched_ions[0]
        if len(matched_ions) == 0:
            return []
        else:
            raise ValueError("不应该存在多个离子")


    def countAtoms(self, formula: str):
        N = len(formula)
        stack = [Counter()]
        i = 0
        while i < N:
            # 当前字符为大写字母的情况
            if formula[i].isupper():
                # 获取原子名称
                idx_start = i
                i += 1
                while i < N and formula[i].islower():  # 找到所有连续的小写字母（即完整原子名）
                    i += 1
                name = formula[idx_start:i]

                # 获取并累加原子数量
                idx_start = i
                while i < N and formula[i].isdigit():
                    i += 1
                stack[-1][name] += int(formula[idx_start:i] or 1)

            # 当前字符为左括号的情况
            elif formula[i] == "(":
                stack.append(Counter())
                i += 1

            elif formula[i] == ")":
                inner = stack.pop()
                i += 1
                idx_start = i
                while i < N and formula[i].isdigit():
                    i += 1
                rate = int(formula[idx_start:i] or 1)
                for name, value in inner.items():
                    stack[-1][name] += value * rate

            else:
                # 跳过非元素、非括号、非数字字符
                logger.error("Invalid character encountered in the formula")
                return None

        # 排序并输出
        count = stack[0]
        ans = {}
        for name in sorted(count):
            ans[name] = count[name]
        return ans

    def format_chemical_formula(self, atoms):
        """
        将原子计数字典转换为化学式字符串。如果字典中包含电荷（'e'键），还将添加电荷表示。

        :param atoms: dict
            包含元素符号和对应数量的字典。可包含键 'e'，表示电荷数（正数表示正电荷，负数表示负电荷）。

        :return: str
            格式化后的化学式，可能包含电荷表示。
        """
        formula = ""
        charge = ""

        for element, count in atoms.items():
            if element != 'e':
                formula += f"{element}{count}" if count > 1 else f"{element}"

        if 'e' in atoms:
            charge_magnitude = abs(atoms['e'])
            charge_sign = '+' if atoms['e'] < 0 else '-'
            if charge_magnitude == 1:
                charge = f"{charge_sign}"
            else:
                charge = f"{charge_magnitude}{charge_sign}"

        if charge:
            return f"[{formula}]{charge}"
        else:
            return formula

    def lipidParser(self, formula):

        """
        解析甘油脂分子式
        :param formula: 甘油脂分子式，例如TG (O-50:1) [NL-16:0]
        :return: 甘油分子的元素和原子数目，{'C': 53, 'H': 102, 'O': 5}
        """

        parts = [self.remove_outer_brackets(part) for part in self.split_lipid_string(formula)]

        logger.debug(f"{self.formula} - first split: {parts}")

        deuterium = self.countDeuterium(parts)

        ions = self.countIon(parts)

        nl_chain = self.countNL(parts)

        lipid_formula = [self.detect_lipid_category(part, lipid_class_dict) for part in parts if self.detect_lipid_category(part, lipid_class_dict)[0]]

        lipid_formula = self.split_lipid(lipid_formula)

        lipid_class = lipid_formula[0][0]
        lipid_sub_class = lipid_formula[0][1]
        # 统计脂肪酸链数目
        chain_num = len(lipid_formula[0][3])
        chain_index = [index for index, element in enumerate(lipid_formula[0][3]) if element[0] != '0:0']
        chain_null_index = [index for index, element in enumerate(lipid_formula[0][3]) if element[0] == '0:0']
        logger.debug(f"{self.formula} - Splitted chains: {lipid_formula[0][3]}, chain_index: {chain_index}, chain_null_index: {chain_null_index}")
        methyl_num = 0
        total_chain = []
        plasmalogen_tag = None
        sm_hydroxyl_tag = None
        chain_hydroxyl_num = 0
        for chain in [lipid_formula[0][3][i] for i in chain_index]:
            if chain[1] == 1:
                plasmalogen_tag = "Alkyl"
            elif chain[1] == 2:
                plasmalogen_tag = "Alkenyl"

            if chain[3] in ["m", "d", "t"]:
                sm_hydroxyl_tag = chain[3]
            methyl_num += chain[2]

            # 羟基数目，若该链为鞘脂骨架，则不计入羟基数目
            if chain[3] is None:
                chain_hydroxyl_num += chain[4]
            total_chain.append(chain[0])

        total_formula = {}
        actual_chain_num = None

        # Bone and nums of Fatty acid chains
        if lipid_class == "GL":
            total_formula = self.merge_element_counts(total_formula, glycerol_atoms)
            if lipid_sub_class == "TG":
                actual_chain_num = 3
            elif lipid_sub_class == "DG":
                actual_chain_num = 2
            elif lipid_sub_class == "MG":
                actual_chain_num = 1

        elif lipid_class == "PL":
            if lipid_sub_class in ["PC", "LPC"]:
                total_formula = self.merge_element_counts(total_formula, glycerol_atoms, pc_atoms)
                total_formula = self.delete_element_counts(total_formula, water_atoms)
            if lipid_sub_class in ["PE", "LPE"]:
                total_formula = self.merge_element_counts(total_formula, glycerol_atoms, pe_atoms)
                total_formula = self.delete_element_counts(total_formula, water_atoms)
            if lipid_sub_class in ["PS", "LPS"]:
                total_formula = self.merge_element_counts(total_formula, glycerol_atoms, ps_atoms)
                total_formula = self.delete_element_counts(total_formula, water_atoms)
            if lipid_sub_class in ["PI", "LPI"]:
                total_formula = self.merge_element_counts(total_formula, glycerol_atoms, pi_atoms)
                total_formula = self.delete_element_counts(total_formula, water_atoms)

            if lipid_sub_class in ["PC", "PE", "PS", "PI"]:
                if len(chain_index) == 1 and len(chain_null_index) == 1:
                    actual_chain_num = 1
                else:
                    actual_chain_num = 2
            elif lipid_sub_class in ["LPC", "LPE", "LPS", "LPI"]:
                actual_chain_num = 1

        elif lipid_class == "SP":
            actual_chain_num = 1
            if sm_hydroxyl_tag is not None:
                if len(chain_index) == 2:
                    sph_bone = sph.from_chain(f"{sm_hydroxyl_tag}{total_chain[0]}")
                    total_chain = [total_chain[1]]
                elif len(chain_index) == 1:
                    # 无法区分骨架，假设骨架为16:0
                    sph_bone = sph.from_chain(f"{sm_hydroxyl_tag}16:0")
                    total_chain[0] = self.delete_lipid_chains(total_chain[0], "16:0")
                else:
                    logger.error(f"{self.formula} - 无法识别鞘脂骨架")
                    return None
            else:
                logger.error(f"{self.formula} - 无法识别鞘脂骨架")
                return None

            if lipid_sub_class == "SM":
                total_formula = self.merge_element_counts(total_formula, sph_bone, pc_atoms)
                total_formula = self.delete_element_counts(total_formula, water_atoms)

            if lipid_sub_class in ["Cer", "meCer"]:
                total_formula = self.merge_element_counts(total_formula, sph_bone)

        # Fatty acid chains
        total_fatty_chain = self.calculate_fatty_acid_chain(total_chain, actual_chain_num)
        # 如果只有一个链，直接计算，否则计算后合并
        if len(total_fatty_chain) > 1:
            total_fatty_chain = self.merge_element_counts(total_fatty_chain)
        else:
            total_fatty_chain = total_fatty_chain[0]

        total_formula = self.merge_element_counts(total_formula, total_fatty_chain)
        logger.debug(f"{self.formula} - total_fatty_chain: {total_fatty_chain}")

        # 酯化反应，酯键减少氢氧原子数目
        if lipid_class == "GL":
                logger.debug(f"{self.formula} - {lipid_class} - {lipid_sub_class}, water * {actual_chain_num}")
                total_formula = self.delete_element_counts(total_formula, [water_atoms] * actual_chain_num)
        elif lipid_class == "PL":
            if lipid_sub_class in ["PC", "PE", "PS", "PI"]:
                # 双脂肪酸链
                if len(chain_index) == 2 or (len(chain_index) == 1 and len(chain_null_index) == 0):
                    total_formula = self.delete_element_counts(total_formula, [water_atoms] * 2)
                # 溶血磷脂，单脂肪酸链
                elif len(chain_index) == 1 and len(chain_null_index) == 1:
                    total_formula = self.delete_element_counts(total_formula, [water_atoms] * 1)
            # 溶血磷脂，单脂肪酸链
            elif lipid_sub_class in ["LPC", "LPE", "LPS", "LPI"]:
                total_formula = self.delete_element_counts(total_formula, [water_atoms] * 1)
        elif lipid_class in ["SP"]:
            total_formula = self.delete_element_counts(total_formula, [water_atoms] * 1)

        # 加上氘原子数目
        if deuterium:
            total_formula = self.merge_element_counts(total_formula, deuterium)
            total_formula = self.delete_element_counts(total_formula, {'H': deuterium['D']})
            logger.info(f"{self.formula} - 检测到氘原子数目：{deuterium['D']}")

        # 缩醛磷脂标记
        if plasmalogen_tag == "Alkyl":
            total_formula = self.delete_element_counts(total_formula, {'O': 1})
            total_formula = self.merge_element_counts(total_formula, {'H': 2})
        elif plasmalogen_tag == "Alkenyl":
            total_formula = self.delete_element_counts(total_formula, {'O': 1})

        # 加上甲基化标记
        if methyl_num > 0:
            total_formula = self.merge_element_counts(total_formula, {'C': methyl_num, 'H': 2 * methyl_num})

        # 加上羟基标记
        if chain_hydroxyl_num > 0:
            total_formula = self.merge_element_counts(total_formula, {'O': chain_hydroxyl_num})

        # 加上离子
        if len(ions) > 0:
            logger.info(f"{self.formula} - 检测到离子")
            total_formula = self.merge_element_counts(total_formula, ions)


        return total_formula

    def Parser(self, formula):

        """
        解析分子式

        :param formula: 输入的任意分子式、或脂质简写形式
        :return: 返回化学式、精确分子量和各元素原子数目
        """

        from molmass import Formula
        from molmass import FormulaError

        if formula == "" or formula is None:
            return None

        self.formula = formula

        category, abbreviation, _ = self.detect_lipid_category(formula, lipid_class_dict)

        if category in ["FA", "GL", "PL", "SP"]:
            logger.info(f"Lipid belongs to category: {category} with abbreviation: {abbreviation}")
            self.atom_list = self.lipidParser(formula)
            total_formula = self.format_chemical_formula(self.atom_list)
            try:
                return total_formula, Formula(total_formula).monoisotopic_mass, self.atom_list
            except FormulaError:
                logger.error(f"化学式不合规: {total_formula}")
                return formula, None, None

        else:
            try:
                logger.info(f"{formula} - 可能为无机物")
                return formula, self.countAtoms(formula), None
            except FormulaError:
                logger.error(f"化学式不合规: {formula}")
                return formula, None, None



    def get_mass(self, atoms : dict = None):

        if atoms is None:
            if self.atom_list is None:
                return None
            atoms = self.atom_list

        formula = self.format_chemical_formula(atoms)

        from molmass import Formula
        return Formula(formula).monoisotopic_mass



