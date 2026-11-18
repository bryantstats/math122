// Math engine for the "p-Series, Zeta, and Primes" slides.
// Everything here uses real numbers only (no complex arithmetic is needed).

export const EULER_GAMMA = 0.5772156649015329;
export const MERTENS = 0.2614972128476428;

// ---------- Primes ----------
export function sieve(n) {
  const comp = new Uint8Array(n + 1);
  comp[0] = 1; if (n >= 1) comp[1] = 1;
  for (let i = 2; i * i <= n; i++) if (!comp[i]) for (let j = i * i; j <= n; j += i) comp[j] = 1;
  return comp; // comp[k] === 0  <=>  k is prime
}

export function primesUpTo(n) {
  const comp = sieve(n), out = [];
  for (let k = 2; k <= n; k++) if (!comp[k]) out.push(k);
  return out;
}

// ---------- Logarithmic integral ----------
export function li(x) {
  const L = Math.log(x);
  let sum = 0, term = 1;
  for (let k = 1; k < 200; k++) {
    term *= L / k;               // L^k / k!
    const add = term / k;
    sum += add;
    if (add < 1e-17 * sum) break;
  }
  return EULER_GAMMA + Math.log(L) + sum;
}

// ---------- Gamma (Lanczos) ----------
const LG = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
  -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
export function gamma(x) {
  if (x < 0.5) return Math.PI / (Math.sin(Math.PI * x) * gamma(1 - x));
  x -= 1;
  let a = LG[0];
  const t = x + 7.5;
  for (let i = 1; i < 9; i++) a += LG[i] / (x + i);
  return Math.sqrt(2 * Math.PI) * Math.pow(t, x + 0.5) * Math.exp(-t) * a;
}

// ---------- Alternating zeta eta(s), s > 0 (Borwein's acceleration) ----------
const BN = 40;
const BD = (() => {
  const d = [];
  let t = 1, s = 1;          // i = 0 term equals 1
  d.push(s);
  for (let i = 1; i <= BN; i++) {
    t *= (4 * (BN + i - 1) * (BN - i + 1)) / ((2 * i) * (2 * i - 1));
    s += t;
    d.push(s);
  }
  return d;
})();
export function eta(s) {
  let sum = 0;
  for (let k = 0; k < BN; k++) sum += ((k % 2 ? -1 : 1) * (BD[k] - BD[BN])) / Math.pow(k + 1, s);
  return -sum / BD[BN];
}

// ---------- Riemann zeta on the real line ----------
export function zeta(s) {
  if (s === 1) return Infinity;
  if (s > 30) return 1 + Math.pow(2, -s) + Math.pow(3, -s);
  if (s > 0) return eta(s) / (1 - Math.pow(2, 1 - s));
  if (s === 0) return -0.5;
  // functional equation: zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s)
  if (Number.isInteger(s) && s % 2 === 0) return 0;
  return Math.pow(2, s) * Math.pow(Math.PI, s - 1) * Math.sin(Math.PI * s / 2) * gamma(1 - s) * zeta(1 - s);
}

// ---------- Hardy's Z function on the critical line (Riemann–Siegel) ----------
export function theta(t) {
  return (t / 2) * Math.log(t / (2 * Math.PI)) - t / 2 - Math.PI / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3);
}
export function hardyZ(t) {
  const a = Math.sqrt(t / (2 * Math.PI));
  const N = Math.floor(a), p = a - N, th = theta(t);
  let sum = 0;
  for (let n = 1; n <= N; n++) sum += Math.cos(th - t * Math.log(n)) / Math.sqrt(n);
  const c0 = Math.cos(2 * Math.PI * (p * p - p - 1 / 16)) / Math.cos(2 * Math.PI * p);
  return 2 * sum + (N % 2 ? 1 : -1) * Math.pow(a, -0.5) * c0;
}

// Zeros of Z(t) on [tmin, tmax] by scanning + bisection.
export function zerosZ(tmax, tmin = 10, step = 0.02) {
  const out = [];
  let t0 = tmin, z0 = hardyZ(t0);
  for (let t1 = tmin + step; t1 <= tmax; t1 += step) {
    const z1 = hardyZ(t1);
    if (z0 === 0 || z0 * z1 < 0) {
      let lo = t0, hi = t1, flo = z0;
      for (let i = 0; i < 50; i++) {
        const m = (lo + hi) / 2, fm = hardyZ(m);
        if (flo * fm <= 0) hi = m; else { lo = m; flo = fm; }
      }
      out.push((lo + hi) / 2);
    }
    t0 = t1; z0 = z1;
  }
  return out;
}

// ---------- Chebyshev psi and Riemann's explicit formula ----------
export function psiSteps(xmax) {
  // returns [{x, y}] step points of psi(x) = sum_{p^m <= x} ln p
  const comp = sieve(xmax), pts = [{ x: 1, y: 0 }];
  let y = 0;
  for (let n = 2; n <= xmax; n++) {
    let p = 0;
    if (!comp[n]) p = n;
    else {
      // is n a prime power?
      for (let q = 2; q * q <= n; q++) if (n % q === 0) {
        let m = n; while (m % q === 0) m /= q;
        if (m === 1) p = q;
        break;
      }
    }
    if (p) { pts.push({ x: n, y }); y += Math.log(p); pts.push({ x: n, y }); }
  }
  pts.push({ x: xmax, y });
  return pts;
}
export function explicitPsi(x, zeros) {
  let s = x - Math.log(2 * Math.PI) - 0.5 * Math.log(1 - 1 / (x * x));
  const L = Math.log(x), r = Math.sqrt(x);
  for (const g of zeros) s -= (2 * r * (0.5 * Math.cos(g * L) + g * Math.sin(g * L))) / (0.25 + g * g);
  return s;
}

// ---------- Small helpers ----------
export function gcd(a, b) { while (b) [a, b] = [b, a % b]; return a; }
export function linspace(a, b, n) { return Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1)); }
