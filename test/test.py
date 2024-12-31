import unittest
from LipidMol import FormulaParser


class TestInvalidInput(unittest.TestCase):
    def setUp(self):
        self.fp = FormulaParser()

    def test_invalid_input(self):
        """测试无效的化学式输入"""
        result = self.fp.Parser("INVALID_INPUT")
        self.assertIsNone(result[2])

    def test_empty_input(self):
        """测试空输入"""
        result = self.fp.Parser("")
        self.assertIsNone(result)


class TestPL(unittest.TestCase):
    def setUp(self):
        self.fp = FormulaParser()

    def test_pc(self):
        """测试PC，包括Lyso形式"""
        result1 = self.fp.Parser("PC(12:0/20:1(11Z))")
        self.assertEqual(result1[2], {'C': 40, 'H': 78, 'O': 8, 'N': 1, 'P': 1})

        result2 = self.fp.Parser("PC(32:1)")
        self.assertEqual(result2[2], {'C': 40, 'H': 78, 'O': 8, 'N': 1, 'P': 1})

        result3 = self.fp.Parser("PC(19:3(10Z,13Z,16Z)/0:0)")
        self.assertEqual(result3[2], {'C': 27, 'H': 50, 'O': 7, 'N': 1, 'P': 1})

        result4 = self.fp.Parser("LPC(19:3)")
        self.assertEqual(result4[2], {'C': 27, 'H': 50, 'O': 7, 'N': 1, 'P': 1})

        result5 = self.fp.Parser("PC(P-18:1(11Z)/16:0)")
        self.assertEqual(result5[2], {'C': 42, 'H': 82, 'O': 7, 'N': 1, 'P': 1})


class TestGL(unittest.TestCase):
    """测试甘油脂的解析"""
    def setUp(self):
        self.fp = FormulaParser()

    def test_fa_simple(self):
        """测试简单脂肪酸（FA）的解析"""
        result = self.fp.Parser("FA(16:1)")
        # ('C16H30O2', 254.22458020604, {'C': 16, 'H': 30, 'O': 2})
        self.assertEqual(result[2], {'C': 16, 'H': 30, 'O': 2})

    def test_gl(self):
        """测试甘油脂"""
        result1 = self.fp.Parser("TG(16:1/18:1/18:1)")
        self.assertEqual(result1[2], {'C': 55, 'H': 100, 'O': 6})

        result2 = self.fp.Parser("TG(20:5(5Z,8Z,11Z,14Z,17Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z)/20:5(5Z,8Z,11Z,14Z,17Z)) (d5)")
        self.assertEqual(result2[2], {'C': 65, 'H': 89, 'O': 6, 'D': 5})

        result3 = self.fp.Parser("TG(16:1(9Z)/16:1(9Z)/17:1(9Z)) [M+H]")
        self.assertEqual(result3[2], {'C': 52, 'H': 95, 'O': 6, 'e': -1})


if __name__ == "__main__":
    unittest.main()