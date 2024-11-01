# S-AES

## 概述
本仓库包含了S-AES（简化AES）算法的实现，以及相关测试和挑战。S-AES是由圣塔·克拉拉大学的Edward Schaefer教授及其学生开发的，旨在作为教育用途的加密算法，虽然它本身并不是一个安全的加密方案。该算法保留了AES的结构和性质，但使用了更少的参数，使得学习和理解变得更加容易。

## 项目结构
- **S-AES算法实现：**
包含S-AES加密和解密功能的源代码。
实现了密钥加（AddRoundKey）、半字节代替（SubNibbles）、行移位（ShiftRows）、列混淆（MixColumns）等核心变换。
- **GUI解密工具：**
提供了图形用户界面（GUI），支持用户输入16位的明文数据和16位的密钥，输出相应的16位密文。
GUI工具使用Python的Tkinter库或其他图形界面库实现。

## 使用方法

### 环境要求

- Python 3.x
- Tkinter（Python 标准库，自带）

### 运行步骤

1. **克隆或下载该项目**:
   - main.py
   
2. **安装所需的库**（如需要）:
   - Tkinter通常与Python自带，无需额外安装。

3. **运行程序**:
   执行以下命令启动应用：
   python main.py

4. **使用GUI界面**:
   - 输入明文（可以是二进制或字母）。
   - 输入正确长度为10位的密钥。
   - 点击“加密”按钮，以获取明文的二进制密文输出。
   - 点击“解密为二进制”按钮，可以将密文转为二进制并显示。
   - 点击“解密为字母”按钮，将密文解密为可读文本。
   - 使用暴力破解功能，输入明文和密文，尝试找出正确的密钥。


## 注意事项

- 密钥必须是 10 位的二进制数。
- 明文的长度必须是 8 的倍数或完整的ASCII字符。
- 如果输入无效，程序将显示错误消息。

## 许可证

此项目为开源项目，您可以自由使用和修改。希望您在使用过程中能够获得乐趣！
