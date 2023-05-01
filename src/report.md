## *Fast Fourier Transform 알고리즘*

### - FFT란?
FFT(고속 푸리에 변환, Fast Fourier Transform) 알고리즘은 시간 도메인의 신호를 주파수 도메인으로 변환하는 효율적인 알고리즘

`장점` 

**FFT 알고리즘은 DFT(이산 푸리에 변환, Discrete Fourier Transform)의 계산량을 크게 줄여줍니다.** DFT는 주파수 변환을 위해 시간 도메인의 모든 샘플을 사용하는 반면, FFT는 분할 정복(divide and conquer) 방법을 사용하여 계산량을 크게 줄여줍니다.

`알고리즘의 이해` 

FFT 알고리즘의 기본 아이디어는 분할 정복 방법을 사용하여 DFT를 계산하는 것입니다. 주어진 N개의 시간 도메인 샘플을 이용하여 주파수 도메인의 N개의 점을 계산합니다. 이를 위해 입력 신호를 반으로 분할하고, 재귀적으로 작은 크기의 DFT를 계산한 후, 이를 결합하여 전체 DFT를 얻습니다. 분할 정복 알고리즘을 사용하기 때문에 계산 복잡도는 O(N log N)입니다.

`응용분야` 

**주파수 분석, 신호 필터링, 컨볼루션 연산, 스펙트럼 분석, 신호 압축 등에 널리 활용됩니다**. 예를 들어, 음성 신호를 주파수 도메인으로 변환하여 주파수 대역별로 분석하거나, 이미지 처리에서 이미지의 주파수 성분을 추출하는 데 사용될 수 있습니다.

`성능` 

**입력 크기에 따라 달라집니다.** 입력 크기가 크면 계산 복잡도가 증가하여 연산 시간이 증가할 수 있습니다. 그러나 FFT 알고리즘은 DFT에 비해 계산량이 효율적이므로 대부분의 경우 실시간 처리에 적합합니다. 

`알고리즘 구현 전 해석`

FFT 알고리즘은 빠른 계산을 위해 분할정복 방식을 사용합니다. 입력 신호를 N개의 작은 신호로 분할하고, 이 작은 신호에 대해 재귀적으로 FFT 알고리즘을 적용합니다. 이 과정은 입력 신호의 길이 N이 2의 거듭제곱일 때에 가장 효과적입니다.

`알고리즘 구현`

```import java.util.Arrays;

public class FFT {
// x 배열과 y 배열의 크기는 동일하다고 가정합니다.
public static Complex[] fft(Complex[] x) {
int n = x.length;

        // base case
        if (n == 1) {
            return new Complex[] { x[0] };
        }

        // radix 2 Cooley-Tukey FFT
        if (n % 2 != 0) {
            throw new IllegalArgumentException("n must be a power of 2");
        }

        // fft of even terms
        Complex[] even = new Complex[n/2];
        for (int k = 0; k < n/2; k++) {
            even[k] = x[2*k];
        }
        Complex[] q = fft(even);

        // fft of odd terms
        Complex[] odd  = even;  // reuse the array
        for (int k = 0; k < n/2; k++) {
            odd[k] = x[2*k + 1];
        }
        Complex[] r = fft(odd);

        // combine
        Complex[] y = new Complex[n];
        for (int k = 0; k < n/2; k++) {
            double kth = -2 * k * Math.PI / n;
            Complex wk = new Complex(Math.cos(kth), Math.sin(kth));
            y[k]       = q[k].plus(wk.times(r[k]));
            y[k + n/2] = q[k].minus(wk.times(r[k]));
        }
        return y;
    }

    public static void main(String[] args) {
        Complex[] x = {new Complex(1, 0), new Complex(2, 0), new Complex(3, 0), new Complex(4, 0)};
        Complex[] y = fft(x);

        System.out.println(Arrays.toString(y));
    }
}

class Complex {
private final double re;   // the real part
private final double im;   // the imaginary part

    // create a new object with the given real and imaginary parts
    public Complex(double real, double imag) {
        re = real;
        im = imag;
    }

    // return a string representation of the invoking Complex object
    public String toString() {
        if (im == 0) return re + "";
        if (re == 0) return im + "i";
        if (im <  0) return re + " - " + (-im) + "i";
        return re + " + " + im + "i";
    }

    // return abs/modulus/magnitude
    public double abs() {
        return Math.hypot(re, im);
    }

    // return angle/phase/argument, normalized to be between -pi and pi
    public double phase() {
        return Math.atan2(im, re);
    }

    // return a new Complex object whose value is (this + b)
    public Complex plus(Complex b) {
        Complex a = this;             // invoking object
        double real = a.re + b.re;
        double imag = a.im + b.im;
        return new Complex(real, imag);
    }

    // return a new Complex object whose value is (this - b)
    public Complex minus(Complex b) {
        Complex a = this;
        double real = a.re - b.re;
        double imag = a.im - b.im;
        return new Complex(real, imag);
    }
 ```
`특정 응용에 적용`

저는  복소수 배열에 대한 FFT알고리즘을 사용하였습니다.
복소수 배열 x에 대한 FFT를 수행하고, 결과를 출력합니다. Complex 클래스는 복소수를 표현하기 위한 클래스로, 덧셈, 뺄셈, 곱셈 연산이 구현되어 있습니다. fft 메서드는 재귀적으로 FFT 알고리즘을 구현하며, 분할 정복을 이용하여 입력 배열을 반으로 나누고 다시 재귀적으로 FFT를 수행합니다

````
import java.util.Arrays;

public class FFT {


    // 복소수 클래스
    static class Complex {
        double real, imag;

        public Complex(double real, double imag) {
            this.real = real;
            this.imag = imag;
        }

        // 덧셈 연산
        public Complex add(Complex c) {
            return new Complex(this.real + c.real, this.imag + c.imag);
        }

        // 뺄셈 연산
        public Complex subtract(Complex c) {
            return new Complex(this.real - c.real, this.imag - c.imag);
        }

        // 곱셈 연산
        public Complex multiply(Complex c) {
            double real = this.real * c.real - this.imag * c.imag;
            double imag = this.real * c.imag + this.imag * c.real;
            return new Complex(real, imag);
        }
    }

    // FFT 알고리즘 구현
    public static Complex[] fft(Complex[] x) {
        int N = x.length;

        // 재귀 종료 조건
        if (N == 1) {
            return new Complex[] { x[0] };
        }

        // 홀수 인덱스와 짝수 인덱스 분리
        Complex[] even = new Complex[N / 2];
        Complex[] odd = new Complex[N / 2];
        for (int i = 0; i < N / 2; i++) {
            even[i] = x[2 * i];
            odd[i] = x[2 * i + 1];
        }

        // 분할 정복을 이용한 FFT
        Complex[] yEven = fft(even);
        Complex[] yOdd = fft(odd);

        // FFT 결과 계산
        Complex[] y = new Complex[N];
        for (int k = 0; k < N / 2; k++) {
            double angle = -2 * k * Math.PI / N;
            Complex wk = new Complex(Math.cos(angle), Math.sin(angle));
            y[k] = yEven[k].add(wk.multiply(yOdd[k]));
            y[k + N / 2] = yEven[k].subtract(wk.multiply(yOdd[k]));
        }

        return y;
    }

    public static void main(String[] args) {
        Complex[] x = {
            new Complex(1, 0),
            new Complex(2, 0),
            new Complex(3, 0),
            new Complex(4, 0)
        };

        Complex[] y = fft(x);

        System.out.println("FFT 결과:");
        System.out.println(Arrays.toString(y));
    }
}
````

`성능 평가에 대한 그래프`
네이버 검색을 통해 성능 비교 그래프를 하이퍼링크를 통해 넣어보았습니다.
[네이버](https://www.google.com/url?sa=i&url=https%3A%2F%2Fkr.mathworks.com%2Fhelp%2Fmatlab%2Fref%2Ffft.html&psig=AOvVaw33coNcMGrUEVRQwcvDjk8-&ust=1683018677408000&source=images&cd=vfe&ved=0CBEQjRxqFwoTCJCopKXj0_4CFQAAAAAdAAAAABAE)



