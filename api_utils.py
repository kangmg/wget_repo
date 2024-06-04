import pandas as pd
import urllib
import json
from bs4 import BeautifulSoup
from time import sleep
import random

if __IPYTHON__:
  from tqdm.notebook import tqdm
else:
  from tqdm import tqdm

headers = {'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/125.0.0.0 Safari/537.36'}

# 한글 -> 영어 변환
MOOD_TAGS_CONVERSION = {
    "가벼운"            :   "Light",
    "감상적인"          :   "Sentimental",
    "강렬한"            :   "Intense",
    "강한"              :   "Strong",
    "고양되는"          :   "Uplifting",
    "공격적인"          :   "Aggresive",
    "관능적인"          :   "Sensual",
    "구슬픈"            :   "Melancholy",
    "극적인"            :   "Dramatic",
    "긍정적인"          :   "Positive",
    "꿈꾸는 듯한"       :   "Dreamy",
    "낭만적인"          :   "Romantic",
    "달콤한"            :   "Sweet",
    "따뜻한"            :   "Warm",
    "무서운"            :   "Scary",
    "밝은"              :   "Bright",
    "복잡한"            :   "Complex",
    "부드러운"          :   "Smooth",
    "분위기 있는"       :   "Atmospheric",
    "사나운"            :   "Fierce",
    "사랑스러운"        :   "Lovely",
    "서정적인"          :   "Lyrical",
    "섹시한"            :   "Sexy",
    "슬프고도 아름다운" :   "Bittersweet",
    "슬픈"              :   "Sad",
    "신나는"            :   "Exciting",
    "신비한"            :   "Mysterious",
    "심각한"            :   "Serious",
    "외로운"            :   "Lonely",
    "우울한"            :   "Gloomy",
    "웅장한"            :   "Epic",
    "자신만만한"        :   "Confident",
    "잔잔한"            :   "Calm",
    "재미있는"          :   "Fun",
    "점잖은"            :   "Gentle",
    "정력적인"          :   "Energetic",
    "진지한"            :   "Earnest",
    "질주하는"          :   "Driving",
    "차가운"            :   "Cold",
    "친밀한"            :   "Intimate",
    "쾌활한"            :   "Cheerful",
    "편안한"            :   "Relaxed",
    "행복한"            :   "Happy",
    "향수어린"          :   "Nostalgic",
    "화난"              :   "Angry"
}

# 영어 -> 한글 변환
MOOD_TAGS_CONVERSION_REVERSE = {v: k for k, v in MOOD_TAGS_CONVERSION.items()}

# 감정 태그 리스트
AVAILABLE_MOOD_TAGS_KOR = list(MOOD_TAGS_CONVERSION.keys())
AVAILABLE_MOOD_TAGS_EN  = list(MOOD_TAGS_CONVERSION.values())


def get_API_result(moodTag:str, API_KEY, result_max=500)->dict:
  """
  Last.FM API를 이용한 mood tag 기반 음악 검색 API 함수
  """
  # api 결과 받기
  API_URL = f"https://ws.audioscrobbler.com/2.0/?method=tag.gettoptracks&tag={moodTag}&limit={result_max}&api_key={API_KEY}&format=json"
  with urllib.request.urlopen(API_URL) as response:
    status = response.status
    if status == 200:
      result = response.read().decode()
      json_data = json.loads(result)
      return json_data
    else:
      raise Exception(f"Fail to get API response.\n\nStatus code : {status}")

def get_youtube_url(LAST_FM_URL:str, headers:dict=headers)->str:
  '''
  Description
  -----------
  Last.FM 페이지의 youtube 링크를 추출하는 함수

  Parameters
  ----------
  - LAST_FM_URL(str) : Last.FM 페이지 노래 정보 링크
  - header : request header

  Returns
  -------
  - youtube url(str) : 해당 곡의 youtube 링크

  Usage
  -----
  >>>url = 'https://www.last.fm/music/Coldplay/_/The+Scientist'
  ...
  >>> result = parse_url(url)
  >>> print(result)
  '''
  request = urllib.request.Request(LAST_FM_URL, headers=headers)
  response = urllib.request.urlopen(request)
  status = response.status

  if status == 200:
    # LAST_FM_URL 페이지 결과 bs4로 받기
    html_content = response.read()
    soup = BeautifulSoup(html_content, 'html5lib')
    # youtube 링크 요소 찾기
    element = soup.select_one('#mantle_skin > header > div.header-new-inner > div.header-new-content > div > div > a')
    # youtube 링크 획득
    yt_url = element["href"]
    if yt_url:
      return yt_url
    else:
      raise Exception(f"Fail to get youtube url.\n\nLast.FM url : {LAST_FM_URL}")
  # status가 200이 아니면 에러 반환
  else:
    raise Exception(f"Fail to get web response.\n\nStatus code : {status}")
  
def get_youtube_url(LAST_FM_URL:str, headers:dict=headers)->str:
  '''
  Description
  -----------
  Last.FM 페이지의 youtube 링크를 추출하는 함수

  Parameters
  ----------
  - LAST_FM_URL(str) : Last.FM 페이지 노래 정보 링크
  - header : request header

  Returns
  -------
  - youtube url(str) : 해당 곡의 youtube 링크

  Usage
  -----
  >>>url = 'https://www.last.fm/music/Coldplay/_/The+Scientist'
  ...
  >>> result = parse_url(url)
  >>> print(result)
  '''
  request = urllib.request.Request(LAST_FM_URL, headers=headers)
  response = urllib.request.urlopen(request)
  status = response.status

  if status == 200:
    # LAST_FM_URL 페이지 결과 bs4로 받기
    html_content = response.read()
    soup = BeautifulSoup(html_content, 'html5lib')
    # youtube 링크 요소 찾기
    element = soup.select_one('#mantle_skin > header > div.header-new-inner > div.header-new-content > div > div > a')
    # youtube 링크 획득
    yt_url = element["href"]
    if yt_url:
      return yt_url
    else:
      raise Exception(f"Fail to get youtube url.\n\nLast.FM url : {LAST_FM_URL}")
  # status가 200이 아니면 에러 반환
  else:
    raise Exception(f"Fail to get web response.\n\nStatus code : {status}")
  
def data_formatter(API_result:dict, moodTag:str, save_csv=False, headers:dict=headers)->pd.DataFrame:
  """
  Description
  -----------
  API 결과를 DataFrame으로 변환하는 함수

  Parameters
  ----------
  -  API_result(dict)  : Last FM API 반환 결과
  -  moodTag(str)      : mood tag ( AVAILABLE_MOOD_TAGS_EN 에 포함된 태그 )
  -  save_csv(bool)    : csv 파일로 저장 여부
    - True  : csv 파일로 저장
    - False : DataFrame 반환
  - hedaers : request header

  Returns
  -------
  - DataFrame
  """
  # inner funciton for data parsing
  def parser(single_music:dict, moodTag:str, headers:dict)->dict:
    music_name = single_music["name"]
    duration = single_music["duration"]
    artist = single_music["artist"]["name"]
    lastfm_url = single_music["url"]
    youtube_url = get_youtube_url(lastfm_url,  headers)
    return {"artist":artist, "music_name":music_name, "duration":duration, "lastfm_url":lastfm_url, "youtube_url": youtube_url, "mood_EN":moodTag, "mood_KOR":MOOD_TAGS_CONVERSION_REVERSE[moodTag]}

  music_list = API_result["tracks"]["track"]
  parsed_music_list = []
  failed_list = []
  # 0.7 ~ 2 초 time sleep 후에 youtube link 획득
  print(f"Collecting {moodTag} Data . . .")
  for music in tqdm(music_list):
    try:
      random_sleep_time = random.uniform(0.7, 2)
      sleep(random_sleep_time)
      parsed_music_list.append(parser(music, moodTag, headers))
    except:
      print("\nTime Sleep\n************")
      sleep(30)
      try:
        parsed_music_list.append(parser(music, moodTag, headers))
      except:
        failed_list.append(music)
  # API 결과를 DataFrame으로 변환
  df = pd.DataFrame(parsed_music_list)

  # 반환 df 또는 csv로 저장
  if save_csv:
    df.to_csv(f"./mood_data_csv/{moodTag}.csv", index=False)
    return failed_list
  else:
    return df, failed_list
