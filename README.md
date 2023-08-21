# 설명
Domain knowledge를 활용한 sequential transformation을 이용하여 anomaly detection 수행

raw data를 활용하여 anomaly detection 하기에는 어려움 
따라서 sequential transformation을 통해 data를 변화시켜 이상치 탐지 진행.

기존 sequential transformation 

![image](https://github.com/KimChangHyun-design/S_T-using-domain-knowledge/assets/127087508/e382d4b3-5c33-4fe0-9eac-d1fe631c4163)
![image](https://github.com/KimChangHyun-design/S_T-using-domain-knowledge/assets/127087508/ec84e95c-4555-4c66-a98f-ad9884e2e942)

Domain knowledge를 활용한 sequential transformation

![image](https://github.com/KimChangHyun-design/S_T-using-domain-knowledge/assets/127087508/9b1aba3d-dd7e-4673-8fed-7d96674ed762)

# model = simualtion data(1~6)
각 data 마다 적용된 sequential transformation

![image](https://github.com/KimChangHyun-design/S_T-using-domain-knowledge/assets/127087508/01a67997-9c1b-4839-b873-d711becfb1e9)

#결과

![image](https://github.com/KimChangHyun-design/S_T-using-domain-knowledge/assets/127087508/4cc779f8-0176-4356-9318-0df069bfd25c)


# 기반 논문  
Functional Outlier Detection and Taxonomy by Sequential Transformations
https://www.sciencedirect.com/science/article/pii/S0167947320300517

# 사용방법
devtools::install_github("otsegun/fdaoutlier")
위 명령어로 fdaoutlier 패키지를 설치해주세요.

