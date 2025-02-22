from typing import List, Tuple, Optional

import numpy as np
import scipy.optimize
import scipy.special
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.preprocessing import RobustScaler
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.utils.multiclass import unique_labels
from sklearn.metrics import euclidean_distances
from sklearn.metrics import matthews_corrcoef
from sklearn.svm import SVC, SVR

from modalibpy.stages.base import Annotation


def platt_scaling(y: np.ndarray, y_pred: np.ndarray, a0: float, b0: float, eps: float = 1e-7) -> Tuple[float, float]:

    def loss(a: float, b: float) -> float:
        y_proba = scipy.special.expit(a * y_pred + b)
        y_proba = np.clip(y_proba, eps, 1 - eps)
        return -np.mean(y * np.log(y_proba) + (1 - y) * np.log(1 - y_proba))

    a_bounds = (a0 * 0.5, a0 * 2)
    b_bounds = (b0 - abs(b0), b0 + abs(b0))

    a, b = a0, b0
    for _ in range(10):
        a = scipy.optimize.minimize_scalar(lambda a_: loss(a_, b), bounds=a_bounds).x
        b = scipy.optimize.minimize_scalar(lambda b_: loss(a, b_), bounds=b_bounds).x

    return a, b


def compute_intercept(y: np.ndarray, y_pred: np.ndarray, alpha: np.ndarray, C: np.ndarray, eps: float = 1e-7) -> float:
    sv_mask = np.logical_and(eps < alpha, alpha < C - eps)
    weights = sv_mask
    #weights = alpha
    s_ = np.sum(weights)
    if s_ > 0:
        weights = weights / s_
    return np.sum(weights * (y - y_pred))


def compute_sample_weights(y: np.ndarray) -> np.ndarray:
    y = (y > 0).astype(int)
    sample_weight = np.full(len(y), np.nan)
    n_classes = 2
    for c in range(n_classes):
        sample_weight[y == c] = len(y) / (n_classes * np.sum(y == c))
    assert not np.any(np.isnan(sample_weight))
    return sample_weight


def solve_qp(Q: np.ndarray, p: np.ndarray, A: np.ndarray, b: np.ndarray, lb: np.ndarray, ub: np.ndarray, x0: np.ndarray) -> np.ndarray:
    import qpsolvers
    return qpsolvers.solve_qp(
        Q, p, A=A, b=b, lb=lb, ub=ub, solver='cvxopt', initvals=x0,
        abstol=1e-7, reltol=1e-7
    )


class SVM(BaseEstimator, ClassifierMixin):

    def __init__(
            self,
            C: float = 10000,
            platt_scaling: bool = True,
            random_state: int = 0
    ):
        self.C: float = C
        self.platt_scaling: bool = platt_scaling
        self.random_state: int = random_state

        self.alpha : Optional[np.ndarray] = None
        self.ay: Optional[np.ndarray] = None

        self.intercept: float = 0.
        self.a: float = 1.
        self.b: float = 0.

    def _fit(self, K: np.ndarray, y: np.ndarray, alpha0: Optional[np.ndarray] = None) -> Tuple[np.ndarray, float]:

        # Encode labels as either 1 or -1
        y = 2 * ((y > 0).astype(float) - 0.5)

        n = len(K)

        def loss(alpha: np.ndarray) -> float:
            ay = alpha * y
            obj = np.sum(alpha)
            obj -= 0.5 * np.sum(np.outer(ay, ay) * K)
            return -float(obj)

        def gradient(alpha: np.ndarray) -> np.ndarray:
            ay = alpha * y
            grad_alpha = np.ones_like(alpha)
            grad_alpha -= np.sum(np.outer(y, ay) * K, axis=1)

            return -grad_alpha

        # Add equality constraints
        A = []
        A.append(y)
        A = np.asarray(A)
        constraints = [scipy.optimize.LinearConstraint(A, lb=0, ub=0)]

        # Add bounds on Lagrange multipliers
        ub = self.C * compute_sample_weights(y)
        bounds = scipy.optimize.Bounds(lb=0, ub=ub, keep_feasible=True)

        # Compute initial solution
        if alpha0 is None:
            model = SVC(kernel='precomputed', C=self.C, class_weight='balanced', probability=False, tol=1e-6, random_state=self.random_state)
            model.fit(K, y)
            alpha = np.zeros(n)
            alpha[model.support_] = np.abs(model.dual_coef_.flatten())
            x0 = alpha
            intercept = model.intercept_
        else:
            x0 = alpha0
            intercept = 0

        # Solve optimization problem
        #options = {'disp': False}
        #res = scipy.optimize.minimize(loss, x0, jac=gradient, constraints=constraints, bounds=bounds, method='SLSQP', options=options)
        #alpha = res.x
        alpha = x0
        #Q = K * np.outer(y, y)
        #p = -np.ones(len(y))
        #b = np.zeros(1)
        #alpha = solve_qp(Q, p, A, b, np.zeros(len(ub)), ub, x0)

        # Compute intercept
        ay = alpha * y
        y_pred = np.sum(ay[:, np.newaxis] * K, axis=0)
        intercept = compute_intercept(y, y_pred, alpha, ub)

        return ay, intercept

    def fit(self, K: np.ndarray, y: np.ndarray) -> 'CARSVM':

        # Fit model
        ay, intercept = self._fit(K, y)
        self.ay = ay
        self.alpha = np.abs(ay)
        self.intercept = intercept

        # Platt scaling
        if self.platt_scaling:
            model = SVC(kernel='precomputed', C=self.C, probability=True, tol=1e-6, random_state=self.random_state)
            model.fit(K, y)
            a0, b0 = -model.probA_[0], -model.probB_[0]
            #a, b = platt_scaling(y, y_pred, a0, b0)
            a, b = a0, b0
            self.a = a
            self.b = b

        return self

    def _decision_function(self, K: np.ndarray, ay: np.ndarray, intercept: float) -> np.ndarray:
        y_pred = np.sum(ay[np.newaxis, :] * K, axis=1)
        return y_pred + intercept

    def decision_function(self, K: np.ndarray) -> np.ndarray:
        return self._decision_function(K, self.ay, self.intercept)

    def predict(self, K: np.ndarray) -> np.ndarray:
        return (self.predict_proba(K)[:, 1] > 0).astype(int)

    def predict_proba(self, K: np.ndarray) -> np.ndarray:
        y_proba = scipy.special.expit(self.a * self.decision_function(K) + self.b)
        return np.asarray([1 - y_proba, y_proba]).T


class CARSVM(BaseEstimator, ClassifierMixin):

    def __init__(
            self,
            C_c: float = 10000,
            C_r: float = 1,
            epsilon: float = 0.2,
            platt_scaling: bool = True,
            balanced: bool = True,
            random_state: int = 0
    ):
        self.C_c: float = C_c
        self.C_r: float = C_r
        self.epsilon: float = epsilon
        self.random_state: int = random_state

        self.alpha_c : Optional[np.ndarray] = None
        self.ay_c: Optional[np.ndarray] = None
        self.intercept_c: float = 0.
        self.a_c: float = 1.
        self.b_c: float = 0.
        self.platt_scaling: bool = platt_scaling
        self.balanced: bool = balanced

        self.scaler_r: RobustScaler = RobustScaler()

    def _fit(
            self,
            K: np.ndarray,
            y_c: np.ndarray,
            y_r: np.ndarray,
            x0: Optional[Tuple[np.ndarray, np.ndarray, np.ndarray]] = None
    ) -> Tuple[np.ndarray, float]:

        n = len(K)

        # Check missing values
        mask_r = ~np.isnan(y_r)
        y_r[~mask_r] = 0

        # Add equality constraints
        A = []
        A.append(np.concatenate((y_c, np.zeros(2 * n)), axis=0))
        A.append(np.concatenate((np.zeros(n), np.ones(n), -np.ones(n)), axis=0))
        A = np.asarray(A)
        constraints = [scipy.optimize.LinearConstraint(A, lb=0, ub=0)]

        # Add bounds on Lagrange multipliers
        sample_weight = compute_sample_weights(y_c)
        ub_c = self.C_c * sample_weight
        ub_rp = self.C_r * sample_weight * mask_r
        ub_rn = self.C_r * sample_weight * mask_r
        ub = np.concatenate((ub_c, ub_rp, ub_rn), axis=0)
        lb = np.zeros_like(ub)
        bounds = scipy.optimize.Bounds(lb=0, ub=ub, keep_feasible=True)

        # Compute initial solution
        if x0 is None:
            model_c = SVC(kernel='precomputed', C=self.C_c, class_weight='balanced', probability=False, tol=1e-6, random_state=self.random_state)
            model_c.fit(K, (y_c > 0).astype(int))
            alpha_c = np.zeros(n)
            alpha_c[model_c.support_] = np.abs(model_c.dual_coef_.flatten())
            alpha_r = np.zeros(n)
            alpha_r[mask_r] = self.fit_svr(K[mask_r, :][:, mask_r], y_r[mask_r])
            alpha_rp = np.copy(alpha_r)
            alpha_rn = np.zeros_like(alpha_rp)
            x0 = np.concatenate((alpha_c, alpha_rp, alpha_rn), axis=0)
        else:
            x0 = np.concatenate(x0, axis=0)

        Q = np.zeros((3 * n, 3 * n))
        Q[:n, :n] = K * np.outer(y_c, y_c)
        Q[n:2*n, n:2*n] = K
        Q[2*n:3*n, 2*n:3*n] = K
        Q[n:2*n, 2*n:3*n] = -K
        Q[2*n:3*n, n:2*n] = -K
        Q[:n, n:2*n] = K * y_c[:, np.newaxis]
        Q[n:2*n, :n] = Q[:n, n:2*n].T
        Q[:n, 2*n:3*n] = -K * y_c[:, np.newaxis]
        Q[2*n:3*n, :n] = Q[:n, 2*n:3*n].T

        p = np.concatenate((
            -np.ones(len(y_c)),
            self.epsilon - y_r,
            self.epsilon + y_r
        ), axis=0)
        b = np.zeros(2)
        alpha = solve_qp(Q, p, A, b, lb, ub, x0.flatten())
        alpha_c = alpha[:n]
        alpha_r = alpha[n:2*n] - alpha[2*n:]

        # Compute intercept for the classification task
        ay_c = alpha_c * y_c + alpha_r
        y_pred = np.sum(ay_c[:, np.newaxis] * K, axis=0)
        intercept_c = compute_intercept(y_c, y_pred, alpha_c, ub_c)
        
        # Compute intercept for the regression task
        #y_pred = np.sum(self.alpha_r[:, np.newaxis] * K, axis=0)
        #self.b_r = np.mean(y_r - y_pred) * self.a_c

        return ay_c, alpha_r, intercept_c


    def fit(self, K: np.ndarray, y: np.ndarray) -> 'CARSVM':

        n = len(K)
        y_c = y[:, 0]
        y_r = y[:, 1]

        # Scale regression targets
        mask = ~np.isnan(y_r)
        self.scaler_r.fit(y_r[mask, np.newaxis])
        y_r = self.scaler_r.transform(y_r[:, np.newaxis]).flatten()

        # Encode labels as either 1 or -1
        y_c = 2 * ((y_c > 0).astype(float) - 0.5)

        # Fit model
        ay_c, alpha_r, intercept_c = self._fit(K, y_c, y_r)
        self.intercept_c = intercept_c

        # Store dual solution
        self.alpha_c = np.abs(ay_c)
        self.alpha_r = alpha_r
        self.ay_c = ay_c

        # Platt scaling
        if self.platt_scaling:
            model = SVC(kernel='precomputed', C=self.C_c, probability=True, tol=1e-6, random_state=self.random_state)
            model.fit(K, y_c)
            a, b = -model.probA_[0], -model.probB_[0]
            self.a_c = a
            self.b_c = b

        return self

    def fit_svr(self, K: np.ndarray, y: np.ndarray) -> 'SVR':
        n = len(K)
        alpha_r = np.zeros(n)
        if self.C_r > 0:
            model_r = SVR(kernel='precomputed', epsilon=self.epsilon, C=self.C_r, tol=1e-6)
            model_r.fit(K, y)
            alpha_r[model_r.support_] = model_r.dual_coef_.flatten()
        return alpha_r

    def decision_function(self, K: np.ndarray) -> np.ndarray:
        coef = self.ay_c + self.alpha_r
        y_pred = np.sum(coef[np.newaxis, :] * K, axis=1)
        return y_pred + self.intercept_c

    def predict(self, K: np.ndarray) -> np.ndarray:
        return (self.predict_proba(K)[:, 1] > 0).astype(int)

    def predict_proba(self, K: np.ndarray) -> np.ndarray:
        y_proba = scipy.special.expit(self.a_c * self.decision_function(K) + self.b_c)
        return np.asarray([1 - y_proba, y_proba]).T


class MultiStageSVM(BaseEstimator, ClassifierMixin):

    def __init__(
            self,
            C: float = 1,
            C_aux: float = 1,
            platt_scaling: bool = True,
            balanced: bool = True,
            random_state: int = 0
    ):
        self.stages: Tuple[..., int] = (0, 1, 2, 3, 4)
        self.C: np.ndarray = np.full(len(self.stages), C_aux)
        self.C[0] = C
        self.models = []
        self.alpha = None
        self.ay = None
        self.intercept = np.zeros(self.n_stages)
        self.a = np.ones(self.n_stages)
        self.b = np.zeros(self.n_stages)
        self.platt_scaling: bool = platt_scaling
        self.balanced: bool = balanced
        self.random_state: int = random_state

    def encode_labels(self, annotations: List[Annotation]) -> np.ndarray:
        y = []
        for annot in annotations:
            grade = float(annot)
            if not annot.is_cancer:
                y.append([-1] * self.n_stages)
            else:
                if np.isnan(grade):
                    y.append([1] + [0] * (self.n_stages - 1))
                else:
                    found = False
                    for k, stage in enumerate(self.stages):
                        if stage <= grade < stage + 1:
                            y.append([1] * (k + 1) + [-1] * (self.n_stages - k - 1))
                            found = True
                            break
                    if not found:
                        raise ValueError(f'Unknown cancer grade "{grade}"')
            assert len(y[-1]) == self.n_stages
        return np.asarray(y, dtype=float)

    def fit_(self, K: np.ndarray, annotations: List[Annotation]) -> 'MultiStageSVM':

        stages = set()
        for annot in annotations:
            if annot.is_cancer:
                stage = float(annot)
                if not np.isnan(stage):
                    stages.add(int(stage))
        self.stages = tuple(sorted(list(stages)))

        n = len(K)
        self.y_signed = self.encode_labels(annotations)
        
        n_classes = self.y_signed.shape[1]

        self.ay = np.zeros((n, n_classes), dtype=float)
        self.is_sv = np.zeros((n, n_classes), dtype=bool)

        self.models = []
        for c in range(n_classes):
            mask = (self.y_signed[:, c] != 0)
            y = (self.y_signed[:, c] > 0).astype(int)
            model = SVC(kernel='precomputed', C=self.C[c], class_weight='balanced', probability=self.platt_scaling, tol=1e-6, random_state=self.random_state)
            if np.any(y[mask]) and np.any(1 - y[mask]):
                model.fit(K[mask, :][:, mask], y[mask])

                v = np.zeros(int(np.sum(mask)))
                v[model.support_] = model.dual_coef_.flatten()
                self.ay[mask, c] += v
                self.is_sv[mask, c] = True

            self.models.append(model)

        self.alpha = np.abs(self.ay)
        self.ay = np.sum(self.ay, axis=1)

        self.refine_svm_solution(K)

        # Infer Platt smoothing parameters
        self.b = np.zeros(self.n_stages)
        self.a = np.ones(self.n_stages)
        if self.platt_scaling:
            for c in range(n_classes):
                if hasattr(self.models[c], 'probA_'):
                    self.a[c] = -self.models[c].probA_
                if hasattr(self.models[c], 'probB_'):
                    self.b[c] = -self.models[c].probB_

        return self

    def fit(self, K: np.ndarray, annotations: List[Annotation]) -> 'MultiStageSVM':
        self.fit_(K, annotations)
        return self

    def refine_svm_solution(self, K: np.ndarray) -> None:

        # Define sample-specific regularization parameters
        lb = 0
        ub = np.full((len(K), self.n_stages), float(self.C[0]))
        for c in range(self.n_stages):
            ub[:, c] = float(self.C[c])
        if self.balanced:
            for c in range(self.n_stages):
                mask = (self.y_signed[:, c] == 0)
                n_known = len(mask) - np.sum(mask)
                ub[mask, c] *= 0
                mask = (self.y_signed[:, c] == 1)
                ub[mask, c] *= n_known / (2 * np.sum(mask))
                mask = (self.y_signed[:, c] == -1)
                ub[mask, c] *= n_known / (2 * np.sum(mask))
        lb = np.zeros_like(ub)

        n = len(K)
        n_stages = self.n_stages

        # Add equality constraints
        A = []
        for c in range(self.y_signed.shape[1]):
            mask = np.zeros(self.y_signed.shape, dtype=bool)
            mask[:, c] = 1
            A.append((self.y_signed * mask).flatten())
        A = np.asarray(A)

        # Solve QP
        Q = np.zeros((self.n_stages * n, self.n_stages * n))
        for c in range(self.n_stages):
            Q[c*n:(c+1)*n, c*n:(c+1)*n] = np.outer(self.y_signed[:, c], self.y_signed[:, c]) * K
        p = -np.ones(self.n_stages * n)
        alpha = solve_qp(Q, p, A, np.zeros(self.n_stages), lb.flatten(), ub.flatten(), self.alpha.flatten())
        alpha = alpha.reshape(n, self.n_stages)
        self.alpha = alpha

        n_classes = self.y_signed.shape[1]
        for c in range(n_classes):
            mask = (self.y_signed[:, c] == 0)
            self.alpha[mask, c] = 0

        assert np.all(lb <= alpha + 1e-5)
        assert np.all(ub >= alpha - 1e-5)

        self.ay = np.sum(self.alpha * self.y_signed, axis=1)

        # Compute intercept
        y_pred = np.sum(self.ay[np.newaxis, :] * K, axis=1)
        for c in range(self.y_signed.shape[1]):
            self.intercept[c] = compute_intercept(self.y_signed[:, c], y_pred, self.alpha[:, c], ub[:, c])

    def decision_function(self, K: np.ndarray) -> np.ndarray:
        y_pred = np.sum(self.ay[np.newaxis, :] * K, axis=1)
        return y_pred + self.intercept[0]

    def predict(self, K: np.ndarray) -> np.ndarray:
        return (self.predict_proba(K)[:, 1] > 0.5).astype(int)

    def predict_proba(self, K: np.ndarray) -> np.ndarray:
        y_proba = scipy.special.expit(self.a[0] * self.decision_function(K) + self.b[0])
        return np.asarray([1 - y_proba, y_proba]).T

    @property
    def n_stages(self) -> int:
        return max(1, len(self.stages))
