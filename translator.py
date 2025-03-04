import openai
import requests



def translator(text):
    """
    使用 OpenAI 或 DeepL API 进行翻译。
    - [0]: 使用 DeepL API
    - [1]: 使用 OpenAI API
    """
    model=['deepl','openai'][1]
    
    if model == 'deepl':
        return deepl_translate(text)
    elif model == 'openai':
        return openai_translate(text)

def openai_translate(text):
    """使用 OpenAI API 进行翻译"""
    openai.api_key = ''  # 你的 OpenAI API Key
    openai.api_base = ''  # 你的 OpenAI API Base URL
    
    try:
        response = openai.ChatCompletion.create(
            model="gpt-3.5-turbo-ca",
            messages=[
                {"role": "system", "content": "你是专业学术翻译，请准确翻译下面的英文文本为中文，保留专有名词、技术术语和公式。"},
                {"role": "user", "content": text}
            ],
            temperature=0.1,
            max_tokens=2048,
            timeout=30
        )
        return response.choices[0].message['content'].strip()
    except Exception as e:
        print(f"OpenAI 翻译错误: {e}")
        return "Translation Error"

def deepl_translate(text):
    """使用 DeepL API 进行翻译"""
    deepl_api_key = ''  # 你的 DeepL API Key
    url = "https://api-free.deepl.com/v2/translate"
    params = {
        "auth_key": deepl_api_key,
        "text": text,
        "target_lang": "ZH"
    }
    
    try:
        response = requests.post(url, data=params)
        response_json = response.json()
        return response_json["translations"][0]["text"]
    except Exception as e:
        print(f"DeepL 翻译错误: {e}")
        return "Translation Error"
