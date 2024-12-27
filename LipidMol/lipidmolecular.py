from collections import Counter, defaultdict
import re
from .lipid_class import lipid_class_dict
from .common_molecular import *



class FormulaParser:

    def __init__(self):
        self.atom_list = None
        self.mass = 0

    def split_lipid_string(self, s):
        # # 检查是否存在以括号结尾后跟有标记的模式，如 "(d5)"
        # match = re.match(r"^(.+?)\s+\((d\d+)\)$", input_string)
        # if match:
        #     # 如果找到匹配，分割为主要部分和标记
        #     return [match.group(1), f"({match.group(2)})"]
        # else:
        #     # 如果没有找到特定模式，返回原字符串作为单个元素的列表
        #     return [input_string]

        # 匹配括号内的内容，包括嵌套的括号，同时保留不在括号内的部分
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
                        if not stack:  # 当栈为空时，添加括号之间的内容
                            parts.append(s[last_index:open_index])
                            parts.append(s[open_index:i + 1])
                            last_index = i + 1
                    elif open_char == '[' and char == ']':
                        if not stack:
                            parts.append(s[last_index:open_index])
                            parts.append(s[open_index:i + 1])
                            last_index = i + 1
                    # 根据需要可以继续处理其他类型的括号，例如大括号等

        # 添加最后一个闭合括号后面的内容
        if last_index < len(s):
            parts.append(s[last_index:])

        # 过滤空字符串和去除多余空格
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
                # 如果参数是字典，直接将其元素和计数加到结果中
                for key, value in arg.items():
                    result[key] += value
            elif isinstance(arg, list):
                # 如果参数是列表，假设列表中包含的都是字典
                for dictionary in arg:
                    for key, value in dictionary.items():
                        result[key] += value
            else:
                raise ValueError(
                    "Unsupported argument type. Each argument must be a dictionary or a list of dictionaries.")

        # 返回一个普通字典格式的结果
        return dict(result)

    def delete_element_counts(self, main_dict, *args):
        # 遍历传入的额外参数
        for arg in args:
            if isinstance(arg, dict):
                # 如果参数是字典，从主字典中减去其元素的计数
                for key, value in arg.items():
                    main_dict[key] = main_dict.get(key, 0) - value
                    # 如果减后的元素值小于等于0，则可以选择删除该键或保留为0
                    if main_dict[key] <= 0:
                        del main_dict[key]
            elif isinstance(arg, list):
                # 如果参数是列表，假设列表中包含的都是字典
                for dictionary in arg:
                    for key, value in dictionary.items():
                        main_dict[key] = main_dict.get(key, 0) - value
                        # 同样的删除或保留逻辑
                        if main_dict[key] <= 0:
                            del main_dict[key]
            else:
                raise ValueError(
                    "Unsupported argument type. Each argument must be a dictionary or a list of dictionaries.")

        return main_dict

    def remove_outer_brackets(self, s):
        # Check for both English and Chinese brackets
        if (s.startswith('(') and s.endswith(')')) or (s.startswith('（') and s.endswith('）')):
            return s[1:-1]
        return s

    def detect_lipid_category(self, lipid_string, lipid_class):
        # Use a regular expression to extract the lipid abbreviation more flexibly
        # It extracts the first word or group of letters before any space or punctuation
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

        return [rc1, rc2, rc3]

    def calculate_fatty_acid_chain(self, chain_list: str | list, lipid_sub_class=None) -> list:

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

        # 如果是TG或DG，需要补充氧原子数目
        if lipid_sub_class == "TG" and len(chain_list) == 1:
            output_chain_list[0]['O'] += 2 * (3 - 1)
        elif lipid_sub_class == "DG" and len(chain_list) == 1:
            output_chain_list[0]['O'] += 2 * (2 - 1)

        return output_chain_list

    def add_lipid_chains(self, chain1, chain2):
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

    def split_lipid(self,formula_list: tuple | list):

        """

        :param formula_list: [('GL', 'TG', 'formula')]
        :return:
        """

        if isinstance(formula_list, tuple):
            formula_list = [formula_list]

        res = []

        for formula in formula_list:
            bone = formula[1]
            formula_str = formula[2]
            # print(f"formula_str: {formula_str}")

            if formula_str.startswith(bone):
                chains = self.remove_outer_brackets(formula_str[len(bone):])
                # print(f"chains: {chains}")
                chain_list = re.split(r'[/_]', chains)
                chain_list = [self.check_chain(chain) for chain in chain_list]
                res.append((formula[0], formula[1], formula[2], chain_list))
            else:
                print(f"ERROR: {formula_str} is not a valid lipid formula")

        return res

    ## 计算氘原子数目
    def countDeuterium(self, formula_list: str | list):
        # Ensure input is a list
        if isinstance(formula_list, str):
            formula_list = [formula_list]

        # List to hold all matched numbers
        deuterium_counts = []

        # Search for 'd' followed by digits in each formula
        for formula in formula_list:
            matches = re.finditer(r'\bd(\d+)\b', formula)
            for match in matches:
                number = match.group(1)
                deuterium_counts.append(int(number))  # Convert string to int and add to list

        return {'D': sum(deuterium_counts)} if deuterium_counts else {}

    def countNL(self, formula_list: str | list):
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
                print("Invalid character encountered in the formula")
                return None

        # 排序并输出
        count = stack[0]
        ans = {}
        for name in sorted(count):
            ans[name] = count[name]
        return ans

    def countLipid(self, formula: str):
        """
        解析脂质分子式

        使用简写的脂质分子式，例如：PC(16:0/18:1)识别脂质

        支持的脂质类别：
            - GL： 甘油脂
                - MG
                - DG
                - TG
            - PL： 磷脂
                - PC
                - PE
                - PS
                - PI
                - PA
            - SP： 鞘脂
                - SM
                - Cer
                - HexCer
            - ST： 固醇
                - ChE 固醇脂

        :param formula:

        :return:
        """
        gl_re = r"^(MG|DG|TG)\b"
        if re.match(gl_re, formula):
            self.atom_list = self.glParser(formula)

    def glParser(self,formula):

        parts = [self.remove_outer_brackets(part) for part in self.split_lipid_string(formula)]

        print(f"parts: {parts}")

        deuterium = self.countDeuterium(parts)

        nl_chain = self.countNL(parts)

        print(f"nl_chain: {nl_chain}")

        lipid_formula = [self.detect_lipid_category(part, lipid_class_dict) for part in parts if self.detect_lipid_category(part, lipid_class_dict)[0]]

        lipid_formula = self.split_lipid(lipid_formula)

        lipid_class = lipid_formula[0][0]
        lipid_sub_class = lipid_formula[0][1]
        chain_num = len(lipid_formula[0][3])
        methyl_num = 0
        total_chain = []
        plasmalogen_tag = None
        for i in range(chain_num):
            chain = lipid_formula[0][3][i]
            if chain[1] == 1:
                plasmalogen_tag = "Alkyl"
            elif chain[1] == 2:
                plasmalogen_tag = "Alkenyl"
            methyl_num += chain[2]
            total_chain.append(chain[0])

        total_formula = {}

        # Bone
        if lipid_class == "GL":
            total_formula = self.merge_element_counts(total_formula, glycerol_atoms)

        # Fatty chain
        total_fatty_chain = self.calculate_fatty_acid_chain(total_chain, lipid_sub_class)
        if len(total_fatty_chain) > 1:
            total_fatty_chain = self.merge_element_counts(total_fatty_chain)
        else:
            total_fatty_chain = total_fatty_chain[0]

        total_formula = self.merge_element_counts(total_formula, total_fatty_chain)

        # 酯化反应，酯键减少氢氧原子数目
        if lipid_sub_class == "TG":
            total_formula = self.delete_element_counts(total_formula, [water_atoms] * 3)
        elif lipid_sub_class == "DG":
            total_formula = self.delete_element_counts(total_formula, [water_atoms] * 2)
        elif lipid_sub_class == "MG":
            total_formula = self.delete_element_counts(total_formula, [water_atoms] * 1)

        # 加上氘原子数目
        if deuterium:
            total_formula = self.merge_element_counts(total_formula, deuterium)
            # print(['H'] * deuterium['D'])
            total_formula = self.delete_element_counts(total_formula, {'H': deuterium['D']})
            print(f"检测到氘原子数目：{deuterium['D']}")

        # 缩醛磷脂标记
        if plasmalogen_tag == "Alkyl":
            total_formula = self.delete_element_counts(total_formula, {'O': 1})
            total_formula = self.merge_element_counts(total_formula, {'H': 2})
        elif plasmalogen_tag == "Alkenyl":
            total_formula = self.delete_element_counts(total_formula, {'O': 1})

        # 加上甲基化标记
        if methyl_num > 0:
            total_formula = self.merge_element_counts(total_formula, {'C': methyl_num, 'H': 2 * methyl_num})

        total_mass = self.get_mass(total_formula)

        # print(f"Total mass: {total_mass}")
        # print(f"total_fatty_chain: {total_fatty_chain}")
        # print(f"Total formula: {total_formula}")

        return total_formula

    def Parser(self, formula, get_mass = False):

        from molmass import Formula
        from molmass import FormulaError

        category, abbreviation, _ = self.detect_lipid_category(formula, lipid_class_dict)

        if category == "GL":
            print(f"Lipid belongs to category: {category} with abbreviation: {abbreviation}")
            self.atom_list = self.glParser(formula)
            total_formula = ''
            for element in self.atom_list.keys():
                count = self.atom_list.get(element, 0)
                total_formula += f"{element}{count}"
            try:
                return total_formula, Formula(total_formula).monoisotopic_mass, self.atom_list
            except FormulaError:
                print(f"Formula is not valid: {total_formula}")
                return formula, None, None

        else:
            try:
                print(f"Formula maybe inorganic")
                return formula, self.countAtoms(formula), None
            except FormulaError:
                print(f"Formula is not valid: {formula}")
                return formula, None, None



    def get_mass(self, atoms : dict = None):
        mass = 0
        if atoms is None:
            if self.atom_list is None:
                return None
            atoms = self.atom_list

        formula = ''
        for element in atoms.keys():
            count = atoms.get(element, 0)
            formula += f"{element}{count}"

        from molmass import Formula
        return Formula(formula).monoisotopic_mass



