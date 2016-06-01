Python Geo Converter
===

파이썬으로 개발 도중 수집한 각종 좌표계를 자유자재로 변환할 일이 있었습니다. 그러나, 다음 API의 경우 일일 3만회의 제한이
있어서 빅데이터 처리에 굉장히 제약이 많았습니다. 그러던 차에 안드로이드펍에서 좋은 소스를 찾아서 Python으로 변환하였고
다시 공유하게 되었습니다.

소스 하단에 `if __name__ == "__main__"` : 이하로 내용을 적어두었습니다.
모듈로 사용하고자 할때는

```python
import GeoConverter

pt = GeoConverter.GeoPoint(x, y)
output = GeoConverter.convert(GeoConverter.TM, GeoConverter.GEO, pt)
```

이런식으로 사용하시면 됩니다 :)




출처 : http://www.androidpub.com/1043970


**원본소스 상단 :)**

```
/**
 * @author aquilegia
 *
 * The code based on hyosang(http://hyosang.kr/tc/96) and aero's blog ((http://aero.sarang.net/map/analysis.html)
 *	License:  LGPL : http://www.gnu.org/copyleft/lesser.html
 */
 ```
