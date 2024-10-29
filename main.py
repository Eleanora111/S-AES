import tkinter as tk
from tkinter import messagebox
import random

# 轮常数
RCON1 = "10000000"
RCON2 = "00110000"

# S-box
S_Box = [['9', '4', 'A', 'B'],
         ['D', '1', '8', '5'],
         ['6', '2', '0', '3'],
         ['C', 'E', 'F', '7']]

# 逆S-box
inv_S_Box = [['A', '5', '9', 'B'],
             ['1', '7', '8', 'F'],
             ['6', '0', '2', '3'],
             ['C', '4', 'D', 'E']]

# 混淆矩阵和逆混淆矩阵
mixMatrix = [['1', '4'], ['4', '1']]
inv_MixMatrix = [['9', '2'], ['2', '9']]


# 按位异或
def xor(str1, str2):
    if len(str1) != len(str2):
        raise ValueError("Error: Length of str1 and str2 should be equal.")

    result = ""
    for i in range(len(str1)):
        if str1[i] not in "01" or str2[i] not in "01":
            raise ValueError("Error: Both input strings must be binary strings.")

        result += '0' if str1[i] == str2[i] else '1'

    return result


# s盒变换
def sBoxTransform(text, sBox):
    locate = 0
    result = ""
    while locate < len(text):
        row = int(text[locate: locate + 2], 2)
        col = int(text[locate + 2: locate + 4], 2)
        result += sBox[row][col]
        locate += 4

    result = bin(int(result, 16))[2:].zfill(8)
    return result


def G(key, cons, sBox):
    rotatedKey = key[4:] + key[:4]
    result = sBoxTransform(rotatedKey, sBox)
    result = xor(result, cons)
    return result


def keyExpansion(key, cons, sBox):
    rightKey = key[8:]
    leftKey = key[0:8]
    newleftKey = xor(G(rightKey, cons, sBox), leftKey)
    newrightKey = xor(newleftKey, rightKey)
    return newleftKey + newrightKey


def halfByteSubstitute(text, sBox):
    locate = 0
    result = ""
    while locate < len(text):
        row = int(text[locate: locate + 2], 2)
        col = int(text[locate + 2: locate + 4], 2)
        result += sBox[row][col]
        locate += 4

    result = bin(int(result, 16))[2:].zfill(16)
    return result


def toStateMatrix(text):
    return [[text[0:4], text[4:8]], [text[8:12], text[12:16]]]


def reStateMatrix(stateMatrix):
    return stateMatrix[0][0] + stateMatrix[0][1] + stateMatrix[1][0] + stateMatrix[1][1]


def left_shift(text):
    stateMatrix = toStateMatrix(text)
    temp = stateMatrix[0][1]
    stateMatrix[0][1] = stateMatrix[1][1]
    stateMatrix[1][1] = temp
    return reStateMatrix(stateMatrix)


def gf_add(a, b):
    return a ^ b


def gf_multiply(a, b, mod=0b10011):
    p = 0
    while b > 0:
        if b & 1:
            p ^= a
        b >>= 1
        a <<= 1
        if a & 0x10:
            a ^= mod
    return p % 16


def colMix(colMatrix, Text):
    myMatrix = toStateMatrix(Text)
    resMatrix = [["", ""], ["", ""]]
    resMatrix[0][0] = bin(gf_add(gf_multiply(int(colMatrix[0][0], 16), int(myMatrix[0][0], 2)),
                                 gf_multiply(int(colMatrix[0][1], 16), int(myMatrix[0][1], 2))))[2:].zfill(4)

    resMatrix[1][0] = bin(gf_add(gf_multiply(int(colMatrix[0][0], 16), int(myMatrix[1][0], 2)),
                                 gf_multiply(int(colMatrix[0][1], 16), int(myMatrix[1][1], 2))))[2:].zfill(4)

    resMatrix[0][1] = bin(gf_add(gf_multiply(int(colMatrix[1][0], 16), int(myMatrix[0][0], 2)),
                                 gf_multiply(int(colMatrix[1][1], 16), int(myMatrix[0][1], 2))))[2:].zfill(4)

    resMatrix[1][1] = bin(gf_add(gf_multiply(int(colMatrix[1][0], 16), int(myMatrix[1][0], 2)),
                                 gf_multiply(int(colMatrix[1][1], 16), int(myMatrix[1][1], 2))))[2:].zfill(4)

    return reStateMatrix(resMatrix)


def encrypt(plainText, key):
    keylist = [key, keyExpansion(key, RCON1, S_Box), keyExpansion(keyExpansion(key, RCON1, S_Box), RCON2, S_Box)]
    roundKey = keylist[0]
    tempText = xor(plainText, roundKey)
    tempText = halfByteSubstitute(tempText, S_Box)
    tempText = left_shift(tempText)
    tempText = colMix(mixMatrix, tempText)
    roundKey = keylist[1]
    tempText = xor(tempText, roundKey)
    tempText = halfByteSubstitute(tempText, S_Box)
    tempText = left_shift(tempText)
    roundKey = keylist[2]
    tempText = xor(tempText, roundKey)
    return tempText


def decrypt(cipherText, key):
    keylist = [key, keyExpansion(key, RCON1, S_Box), keyExpansion(keyExpansion(key, RCON1, S_Box), RCON2, S_Box)]
    roundKey = keylist[2]
    tempText = xor(cipherText, roundKey)
    tempText = left_shift(tempText)
    tempText = halfByteSubstitute(tempText, inv_S_Box)
    roundKey = keylist[1]
    tempText = xor(tempText, roundKey)
    tempText = colMix(inv_MixMatrix, tempText)
    tempText = left_shift(tempText)
    tempText = halfByteSubstitute(tempText, inv_S_Box)
    roundKey = keylist[0]
    tempText = xor(tempText, roundKey)
    return tempText



# 工具函数：判断输入是否为二进制或字母
def is_binary(input_str):
    return all(c in '01' for c in input_str)

def is_ascii_string(input_str):
    return input_str.isalpha() and all(ord(c) < 128 for c in input_str)

# 字符串转二进制
def string_to_binary(text):
    return ''.join(format(ord(c), '08b') for c in text)

# 二进制转字符串
def binary_to_string(bin_text):
    chars = [chr(int(bin_text[i:i + 8], 2)) for i in range(0, len(bin_text), 8)]
    return ''.join(chars)


# 多重加密
def multi_encrypt(plainText, key):
    # 按位数分割密钥
    if len(key) == 16:
        key_parts = [key]
    elif len(key) == 32:
        key_parts = [key[:16], key[16:]]
    elif len(key) == 48:
        key_parts = [key[:16], key[16:32], key[32:]]
    else:
        return "wrong key"

    # 执行多重加密
    tempText = plainText
    for part in key_parts:
        tempText = encrypt(tempText, part)  # 调用现有的单轮加密

    return tempText


# 中间相遇攻击函数
def meet_in_the_middle_attack(plaintext, ciphertext):
    potential_keys = []

    #构建所有可能的前向和反向加密路径
    for k1 in range(2 ** 16):
        k1_bin = bin(k1)[2:].zfill(16)
        enc_result = encrypt(plaintext, k1_bin)  # 使用k1加密明文
        potential_keys.append((k1_bin, enc_result))

    # Step 2: 使用密文的反向解密路径匹配前向加密路径
    for k2 in range(2 ** 16):
        k2_bin = bin(k2)[2:].zfill(16)
        dec_result = decrypt(ciphertext, k2_bin)  # 使用k2解密密文
        for (k1_candidate, enc_result) in potential_keys:
            if enc_result == dec_result:
                # 找到匹配的密钥对 (k1, k2)
                return k1_candidate + k2_bin  # 返回32位密钥

    return None  # 没有找到匹配的密钥

# 生成随机的16位初始向量
def generate_initial_vector():
    return ''.join(random.choice('01') for _ in range(16))


# 创建 GUI
class SAESApp:
    def __init__(self, master):
        self.master = master
        master.title("S-AES 加密解密")
        master.configure(bg="#F0FFFF")

        # 主框架
        self.frame = tk.Frame(master, bg="#F0FFFF")
        self.frame.pack(padx=20, pady=20)

        # 明文输入
        self.plaintext_label = tk.Label(self.frame, text="明文（16n位二进制或2n个字符）:", bg="#F0FFFF", fg="black",
                                        font=("Arial", 12))
        self.plaintext_label.grid(row=0, column=0, sticky="w")

        self.plaintext_entry = tk.Entry(self.frame, width=30, font=("Arial", 12))
        self.plaintext_entry.grid(row=0, column=1, padx=10, pady=5)

        # 密钥输入
        self.key_label = tk.Label(self.frame, text="密钥（16n位二进制或2n个字符(n=1,2,3)）:", bg="#F0FFFF", fg="black",
                                  font=("Arial", 12))
        self.key_label.grid(row=1, column=0, sticky="w")

        self.key_entry = tk.Entry(self.frame, width=30, font=("Arial", 12))
        self.key_entry.grid(row=1, column=1, padx=10, pady=5)

        # 初始向量IV
        self.iv_label = tk.Label(self.frame, text="初始向量IV（16位二进制）:", bg="#F0FFFF", fg="black",
                                 font=("Arial", 12))
        self.iv_label.grid(row=2, column=0, sticky="w")

        self.iv_entry = tk.Entry(self.frame, width=30, font=("Arial", 12))
        self.iv_entry.grid(row=2, column=1, padx=10, pady=5)

        # 密文输入
        self.ciphertext_label = tk.Label(self.frame, text="密文（16n位二进制或2n个字符）:", bg="#F0FFFF", fg="black",
                                         font=("Arial", 12))
        self.ciphertext_label.grid(row=3, column=0, sticky="w")

        self.ciphertext_entry = tk.Entry(self.frame, width=30, font=("Arial", 12))
        self.ciphertext_entry.grid(row=3, column=1, padx=10, pady=5)

        # 新增按钮
        self.encrypt_bin_button = tk.Button(self.frame, text="加密为二进制", command=self.encrypt_to_binary,
                                            bg="#4CAF50", fg="white", font=("Arial", 12))
        self.encrypt_bin_button.grid(row=4, column=0, padx=10, pady=20)

        self.encrypt_str_button = tk.Button(self.frame, text="加密为字符串", command=self.encrypt_to_string,
                                            bg="#4CAF50", fg="white", font=("Arial", 12))
        self.encrypt_str_button.grid(row=4, column=1, padx=10, pady=20)

        self.decrypt_bin_button = tk.Button(self.frame, text="解密为二进制", command=self.decrypt_to_binary,
                                            bg="#F44336", fg="white", font=("Arial", 12))
        self.decrypt_bin_button.grid(row=5, column=0, padx=10, pady=20)

        self.decrypt_str_button = tk.Button(self.frame, text="解密为字符串", command=self.decrypt_to_string,
                                            bg="#F44336", fg="white", font=("Arial", 12))
        self.decrypt_str_button.grid(row=5, column=1, padx=10, pady=20)

        self.crack_key_button = tk.Button(self.frame, text="破解密钥(32bit)", command=self.crack_key,
                                          bg="#FF9800", fg="white", font=("Arial", 12))
        self.crack_key_button.grid(row=5, column=0, columnspan=2, padx=10, pady=20)

        self.generate_iv_button = tk.Button(self.frame, text="生成初始向量IV", command=self.generate_iv,
                                            bg="#2196F3", fg="white", font=("Arial", 12))
        self.generate_iv_button.grid(row=6, column=0, padx=10, pady=20)

        self.cbc_encrypt_button = tk.Button(self.frame, text="CBC模式加密", command=self.cbc_encrypt,
                                            bg="#4CAF50", fg="white", font=("Arial", 12))
        self.cbc_encrypt_button.grid(row=6, column=1, padx=10, pady=20)

        self.cbc_decrypt_button = tk.Button(self.frame, text="CBC模式解密", command=self.cbc_decrypt,
                                            bg="#F44336", fg="white", font=("Arial", 12))
        self.cbc_decrypt_button.grid(row=6, column=0, columnspan=2, padx=10, pady=20)


    def encrypt_to_binary(self):
        plaintext = self.plaintext_entry.get()
        key = self.key_entry.get()
        if is_ascii_string(plaintext):
            if len(plaintext) % 2 != 0 or len(key) % 2 != 0:
                messagebox.showerror("错误", "明文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            plaintext = string_to_binary(plaintext)
            key = string_to_binary(key)
        if len(plaintext) % 16 != 0 or len(key) % 16 != 0:
            messagebox.showerror("错误", "明文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            return
        ciphertext = ""
        for i in range(0, len(plaintext), 16):
            block = plaintext[i:i + 16]
            ciphertext += multi_encrypt(block, key)
            if ciphertext == "wrong key":
                messagebox.showerror("错误", "密钥长度错误！必须为16位、24位或32位！")
                return
        self.ciphertext_entry.delete(0, tk.END)
        self.ciphertext_entry.insert(0, ciphertext)

    def encrypt_to_string(self):
        plaintext = self.plaintext_entry.get()
        key = self.key_entry.get()
        if is_ascii_string(plaintext):
            if len(plaintext) % 2 != 0 or len(key) % 2 != 0:
                messagebox.showerror("错误", "明文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            plaintext = string_to_binary(plaintext)
            key = string_to_binary(key)
        if len(plaintext) % 16 != 0 or len(key) % 16 != 0:
            messagebox.showerror("错误", "明文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            return
        ciphertext = ""
        for i in range(0, len(plaintext), 16):
            block = plaintext[i:i + 16]
            ciphertext += multi_encrypt(block, key)
            if ciphertext == "wrong key":
                messagebox.showerror("错误", "密钥长度错误！必须为16位、24位或32位！")
                return
        self.ciphertext_entry.delete(0, tk.END)
        self.ciphertext_entry.insert(0, binary_to_string(ciphertext))

    def decrypt_to_binary(self):
        ciphertext = self.ciphertext_entry.get()
        key = self.key_entry.get()
        if is_ascii_string(ciphertext):
            if len(ciphertext) % 2 != 0 or len(key) % 2 != 0:
                messagebox.showerror("错误", "密文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            ciphertext = string_to_binary(ciphertext)
            key = string_to_binary(key)
        if len(ciphertext) % 16 != 0 or len(key) % 16 != 0:
            messagebox.showerror("错误", "密文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            return
        plaintext = ""
        for i in range(0, len(ciphertext), 16):
            block = ciphertext[i:i + 16]
            plaintext += decrypt(block, key)
        self.plaintext_entry.delete(0, tk.END)
        self.plaintext_entry.insert(0, plaintext)

    def decrypt_to_string(self):
        ciphertext = self.ciphertext_entry.get()
        key = self.key_entry.get()
        if is_ascii_string(ciphertext):
            if len(ciphertext) % 2 != 0 or len(key) % 2 != 0:
                messagebox.showerror("错误", "密文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            ciphertext = string_to_binary(ciphertext)
            key = string_to_binary(key)
        if len(ciphertext) % 16 != 0 or len(key) % 16 != 0:
            messagebox.showerror("错误", "密文必须是16位的倍数二进制或2的倍数个字符，密钥必须是16位二进制！")
            return
        plaintext = ""
        for i in range(0, len(ciphertext), 16):
            block = ciphertext[i:i + 16]
            plaintext += decrypt(block, key)
        self.plaintext_entry.delete(0, tk.END)
        self.plaintext_entry.insert(0, binary_to_string(plaintext))

    # 破解密钥功能
    def crack_key(self):
        plaintext = self.plaintext_entry.get()
        ciphertext = self.ciphertext_entry.get()

        # 验证明文和密文格式
        if len(plaintext) != 16 or len(ciphertext) != 16:
            messagebox.showerror("错误", "请输入16位二进制格式的明文和密文！")
            return

        # 尝试中间相遇攻击
        key = meet_in_the_middle_attack(plaintext, ciphertext)

        if key:
            messagebox.showinfo("破解成功", f"找到的密钥为：{key}")
            self.key_entry.delete(0, tk.END)
            self.key_entry.insert(0, key)
        else:
            messagebox.showwarning("破解失败", "未找到匹配的密钥，请检查输入！")

    def generate_iv(self):
        iv = generate_initial_vector()
        self.iv_entry.delete(0, tk.END)
        self.iv_entry.insert(0, iv)

    def cbc_encrypt(self):
        plaintext = self.plaintext_entry.get()
        key = self.key_entry.get()
        iv = self.iv_entry.get()

        if len(iv) != 16 or not is_binary(iv):
            messagebox.showerror("错误", "初始向量IV必须是16位的二进制字符串！")
            return

        if len(plaintext) % 16 != 0 or len(key) % 16 != 0:
            messagebox.showerror("错误", "明文必须是16位的倍数二进制，密钥必须是16位二进制！")
            return

        ciphertext = ""
        previous_block = iv
        for i in range(0, len(plaintext), 16):
            block = plaintext[i:i + 16]
            block = xor(block, previous_block)  # 与前一个块异或
            encrypted_block = encrypt(block, key)  # 加密
            ciphertext += encrypted_block
            previous_block = encrypted_block  # 更新前一个块

        self.ciphertext_entry.delete(0, tk.END)
        self.ciphertext_entry.insert(0, ciphertext)

    def cbc_decrypt(self):
        ciphertext = self.ciphertext_entry.get()
        key = self.key_entry.get()
        iv = self.iv_entry.get()

        if len(iv) != 16 or not is_binary(iv):
            messagebox.showerror("错误", "初始向量IV必须是16位的二进制字符串！")
            return

        if len(ciphertext) % 16 != 0 or len(key) % 16 != 0:
            messagebox.showerror("错误", "密文必须是16位的倍数二进制，密钥必须是16位二进制！")
            return

        plaintext = ""
        previous_block = iv
        for i in range(0, len(ciphertext), 16):
            block = ciphertext[i:i + 16]
            decrypted_block = decrypt(block, key)  # 解密
            plaintext += xor(decrypted_block, previous_block)  # 与前一个块异或
            previous_block = block  # 更新前一个块

        self.plaintext_entry.delete(0, tk.END)
        self.plaintext_entry.insert(0, plaintext)


    # 主程序入口
if __name__ == "__main__":
    root = tk.Tk()
    app = SAESApp(root)
    root.mainloop()