import openai
def translator(text):
    """
    使用 OpenAI API 进行翻译，将英文文本翻译为中文。
    """
    openai.api_key = ''#your key
    openai.api_base = ''#base
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
        translated_text = response.choices[0].message['content'].strip()
    except Exception as e:
        print(f"OpenAI 翻译错误: {e}")
        translated_text = "Translation Error"
    return translated_text