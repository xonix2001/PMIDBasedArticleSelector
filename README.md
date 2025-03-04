# PMIDBasedArticleSelector Usage Guide

## English Instructions

1. **Download or Git Clone:**  
   Download the project or clone it from the repository.
   ```bash
   git clone https://github.com/xonix2001/PMIDBasedArticleSelector.git
   ```

2. **Set Up the Environment:**  
   Install the required packages by running the commands listed in `requirements.txt`.
   ```bash
   pip install -r requirements.txt
   ```

3. **Configure Translator Settings:**  
   Open `translator.py` and enter your OpenAI API key and the base configuration.  
   *If you don't have any key, try to buy one at [ChatAnywhere Shop](https://api.chatanywhere.tech/#/shop).*
   ```python
   # Example in translator.py:
   openai.api_key = '#your key'
   openai.api_base = '#base'
   ```

4. **Prepare PubMed API Account:**  
   Ensure you have a valid PubMed API account.

5. **Enter Your Account Information:**  
   Open `extend.py` and input your account details.
   ```python
   # Example in extend.py:
   email = "#enter here"
   api_key = "#enter here"
   ```

6. **Specify Input and Output Paths:**  
   Set the input path for the PMIDs list (a TXT file) and the desired output path for the resulting Excel file.
   ```python
   # Example usage in extend.py:
   processor_extended.run('#enter here')  # input (.txt) path
   ```

7. **Run the Application:**  
   Execute `extend.py` to start processing.
   ```bash
   python extend.py
   ```

---

## 中文说明

**PMIDBasedArticleSelector 使用指南**

1. **下载或克隆代码：**  
   下载项目或通过 Git 克隆代码仓库。
   ```bash
   git clone https://github.com/xonix2001/PMIDBasedArticleSelector.git
   ```

2. **配置环境：**  
   根据 `requirements.txt` 中的依赖，安装所需的 Python 包。
   ```bash
   pip install -r requirements.txt
   ```

3. **配置翻译器：**  
   打开 `translator.py`，填写你的 OpenAI API 密钥及基础配置。  
   *如果没有 API 密钥，请访问 [ChatAnywhere 商店](https://api.chatanywhere.tech/#/shop) 购买。*
   ```python
   # translator.py 示例:
   openai.api_key = "#your key"
   openai.api_base = "#base"
   ```

4. **准备 PubMed API 账户：**  
   确保你拥有有效的 PubMed API 账户。

5. **填写账户信息：**  
   打开 `extend.py`，输入你的账户相关信息。
   ```python
   # extend.py 示例:
   email = "#enter here"
   api_key = "#enter here"
   ```

6. **设置输入输出路径：**  
   指定包含 PMID 列表的 TXT 文件的输入路径，以及生成结果 Excel 文件的输出路径。
   ```python
   # extend.py 示例:
   processor_extended.run("#enter here")  # 输入文件 (.txt) 路径
   ```

7. **运行程序：**  
   执行 `extend.py` 开始处理。
   ```bash
   python extend.py
   ```
