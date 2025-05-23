import openai
def generate_advanced_pubmed_query(description, date_range=None):
    """
    根据自然语言描述智能生成 PubMed 检索词，并确保日期范围包含且不重复
    :param description: str, 用户提供的自然语言描述
    :param date_range: tuple, 指定日期范围 (start_year, end_year)
    :return: str, 智能生成的 PubMed 检索词
    """
    # 如果有日期范围，生成日期条件
    if date_range and len(date_range) == 2:
        start_year, end_year = date_range
        date_condition = f' AND ("{start_year}"[Date - Publication] : "{end_year}"[Date - Publication])'
    else:
        date_condition = ""

    # GPT 提示词设计
    prompt = f"""
    根据以下描述生成一个 PubMed 检索词：
    描述：{description}
    要求：
    - 自动识别逻辑关系（如 AND, OR, NOT）。
    - 可识别关键词及其优先级，分配合适的 PubMed 检索字段（如 Title、Abstract、MeSH Terms 等）。
    - 如提到排除关系，使用 NOT。
    - 日期范围：{date_range if date_range else "未指定"}。
    - 如果提供了日期范围，请将其转换为 PubMed 格式，例如：
      ("2010"[Date - Publication] : "2024"[Date - Publication])。
    - 不要重复附加日期范围。
    - 返回格式符合 PubMed 检索语法，直接可用。

    示例：
    描述：研究与癫痫相关的神经科学，特别是涉及动作电位，但不包括小鼠实验。
    日期范围：(2010, 2024)
    输出：
    "epilepsy"[Title/Abstract] AND "neuroscience"[Title/Abstract] AND "action potential"[MeSH Terms] NOT "mouse"[Title/Abstract] AND ("2010"[Date - Publication] : "2024"[Date - Publication])
    """
    
    try:
        response = openai.ChatCompletion.create(
            model="gpt-4-turbo-ca",
            messages=[
                {"role": "system", "content": "你是一位熟悉 PubMed 检索语法的助手。"},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7,
            max_tokens=200
        )
    
        query = response['choices'][0]['message']['content'].strip()
        if date_condition not in query:
            query += date_condition
        return query
    except Exception as e:
        print(f"调用 GPT 失败：{e}")
        return None

if __name__ == "__main__":
    # example
    # 设置 OpenAI API Key
    openai.api_key = ""
    openai.api_base = ""
    
    # 用户提供自然语言描述
    description = "研究与癫痫相关的神经科学，尤其是动作电位的作用，但排除关于小鼠的研究。"
    date_range = (2010, 2024)  # 可选日期范围

    # 生成检索词
    query = generate_advanced_pubmed_query(description, date_range)
    print("生成的 PubMed 检索词：")
    print(query)
