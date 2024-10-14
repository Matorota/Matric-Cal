#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

class Fraction {
private:
    // Custom gcd function using the Euclidean algorithm
    int gcd(int a, int b) {
        while (b != 0) {
            int temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }

    // Simplify the fraction by dividing both numerator and denominator by their GCD
    void simplify() {
        int divisor = gcd(numerator, denominator);
        numerator /= divisor;
        denominator /= divisor;
    }

public:
    long long numerator;
    long long denominator;

    Fraction() {
        numerator = denominator = 1;
    }

    // Constructor to initialize the fraction
    Fraction(long long num, long long denom = 1) {
        if (denom == 0) {
            throw std::invalid_argument("Denominator cannot be zero.");
        }
        numerator = num;
        denominator = denom;
        simplify();
    }

    double ToDouble() const {
        return numerator / denominator;
    }

    Fraction abs() {
        return Fraction(numerator < 0 ? numerator * -1 : numerator, denominator < 0 ? denominator * -1 : denominator);
    }

    // Overload arithmetic operators as needed (+, -, *, /)
    Fraction operator+(const Fraction& other) const {
        int num = numerator * other.denominator + other.numerator * denominator;
        int denom = denominator * other.denominator;
        return Fraction(num, denom);
    }

    Fraction& operator+=(const Fraction& other) {
        int num = numerator * other.denominator + other.numerator * denominator;
        int denom = denominator * other.denominator;

        numerator = num;
        denominator = denom;

        simplify();
        return *this;
    }

    Fraction operator-(const Fraction& other) const {
        int num = numerator * other.denominator - other.numerator * denominator;
        int denom = denominator * other.denominator;
        return Fraction(num, denom);
    }

    Fraction& operator-=(const Fraction& other) {
        int num = numerator * other.denominator - other.numerator * denominator;
        int denom = denominator * other.denominator;

        numerator = num;
        denominator = denom;

        simplify();
        return *this;
    }

    Fraction operator*(const Fraction& other) const {
        int num = numerator * other.numerator;
        int denom = denominator * other.denominator;
        return Fraction(num, denom);
    }

    Fraction& operator*=(const Fraction& other) {
        int num = numerator * other.numerator;
        int denom = denominator * other.denominator;

        numerator = num;
        denominator = denom;

        simplify();
        return *this;
    }

    Fraction operator/(const Fraction& other) const {
        if (other.numerator == 0) {
            throw std::invalid_argument("Cannot divide by zero.");
        }
        int num = numerator * other.denominator;
        int denom = denominator * other.numerator;
        return Fraction(num, denom);
    }

    Fraction& operator/=(const Fraction& other) {
        if (other.numerator == 0) {
            throw std::invalid_argument("Cannot divide by zero.");
        }
        int num = numerator * other.denominator;
        int denom = denominator * other.numerator;

        numerator = num;
        denominator = denom;

        simplify();
        return *this;
    }


};

// Overload << operator to display fraction
std::ostream& operator<<(std::ostream& os, const Fraction& f) {
    os << f.numerator;

    if (f.denominator != 1)
        os << "/" << f.denominator;
    return os;
}

std::istream& operator>>(std::istream& is, Fraction& f) {
    long long value;

    is >> value;

    f = Fraction(value);

    return is;
}

static bool operator <(const Fraction lhs, const Fraction rhs) {
    return lhs.ToDouble() < rhs.ToDouble();
}

static bool operator >(const Fraction lhs, const Fraction rhs) {
    return !(lhs < rhs);
}

static bool operator ==(const Fraction lhs, const Fraction rhs) {
    return lhs.ToDouble() == rhs.ToDouble();
}

static bool operator !=(const Fraction lhs, const Fraction rhs) {
    return !(lhs == rhs);
}

static bool operator <=(const Fraction lhs, const Fraction rhs) {
    return (lhs < rhs) || (lhs == rhs);
}

static bool operator >=(const Fraction lhs, const Fraction rhs) {
    return (lhs > rhs) || (lhs == rhs);
}



void printMatrix(const vector<vector<Fraction>> mat) {
    for (const auto row : mat) {
        for (const Fraction val : row) {
            cout << setw(10) << val << "\t";
        }
        cout << endl;
    }
    cout << "----------------------------------\n";
}

void inputMatrix(vector<vector<Fraction>>& matrix, int rows, int cols, string matrixName) {
    cout << "Enter the elements of matrix " << matrixName << " (" << rows << "x" << cols << "):\n";
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            cout << "Enter element [" << i + 1 << "][" << j + 1 << "]: ";
            cin >> matrix[i][j];
        }
    }
}

// Function to create a zero matrix of given dimensions
vector<vector<Fraction>> createZeroMatrix(int rows, int cols) {
    return vector<vector<Fraction>>(rows, vector<Fraction>(cols, 0));
}

void solveMatrixEquation() {
    int n;
    cout << "Enter the number of variables (n): ";
    cin >> n;

    // Initialize augmented matrix (A | B)
    vector<vector<Fraction>> augmented(n, vector<Fraction>(n + 1));

    cout << "\nWe are solving a system of linear equations represented as an augmented matrix (Ax = B).\n";
    cout << "Please input the coefficients of the augmented matrix (n rows, n+1 columns, where the last column is the constant vector B).\n";

    inputMatrix(augmented, n, n + 1, "Augmented Matrix");

    cout << "\nInitial Augmented Matrix:\n";
    printMatrix(augmented);

    // Forward Elimination
    for (int i = 0; i < n; ++i) {
        // Partial Pivoting
        Fraction maxElement = augmented[i][i].abs();
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (augmented[k][i].abs() > maxElement) {
                maxElement = augmented[k][i].abs();
                maxRow = k;
            }
        }

        // Swap maximum row with current row (pivoting)
        if (maxRow != i) {
            swap(augmented[i], augmented[maxRow]);
            cout << "Swapped row " << i + 1 << " with row " << maxRow + 1 << " for pivoting.\n";
            printMatrix(augmented);
        }

        // Make all rows below this one 0 in current column
        for (int k = i + 1; k < n; ++k) {
            Fraction factor = augmented[k][i] / augmented[i][i];
            for (int j = i; j <= n; ++j) {
                (augmented[k][j]) -= factor * augmented[i][j];
            }
            cout << "Eliminated element in row " << k + 1 << ", column " << i + 1 << " using row " << i + 1 << ".\n";
            printMatrix(augmented);
        }
    }

    // Back Substitution
    vector<Fraction> x(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        if (augmented[i][i] == 0) {
            cout << "No unique solution exists.\n";
            return;
        }
        x[i] = augmented[i][n];
        for (int j = i + 1; j < n; ++j) {
            x[i] -= augmented[i][j] * x[j];
        }
        x[i] /= augmented[i][i];
    }

    // Display the solution
    cout << "Solution:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << fixed << setprecision(4) << x[i] << endl;
    }
}

// 2. Matrix Operations (Addition, Subtraction, Multiplication, Division)
vector<vector<Fraction>> addMatrices(const vector<vector<Fraction>>& A, const vector<vector<Fraction>>& B, int r1, int c1, int r2, int c2) {
    int rows = max(r1, r2);
    int cols = max(c1, c2);
    vector<vector<Fraction>> result = createZeroMatrix(rows, cols);

    // Add elements of A
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c1; ++j) {
            result[i][j] += A[i][j];
        }
    }

    // Add elements of B
    for (int i = 0; i < r2; ++i) {
        for (int j = 0; j < c2; ++j) {
            result[i][j] += B[i][j];
        }
    }

    return result;
}

vector<vector<Fraction>> subtractMatrices(const vector<vector<Fraction>>& A, const vector<vector<Fraction>>& B, int r1, int c1, int r2, int c2) {
    int rows = max(r1, r2);
    int cols = max(c1, c2);
    vector<vector<Fraction>> result = createZeroMatrix(rows, cols);

    // Subtract elements of B from A
    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c1; ++j) {
            result[i][j] += A[i][j];
        }
    }

    for (int i = 0; i < r2; ++i) {
        for (int j = 0; j < c2; ++j) {
            result[i][j] -= B[i][j];
        }
    }

    return result;
}

vector<vector<Fraction>> divideMatrices(const vector<vector<Fraction>>& A, const vector<vector<Fraction>>& B, int r1, int c1, int r2, int c2) {
    int rows = max(r1, r2);
    int cols = max(c1, c2);
    vector<vector<Fraction>> result = createZeroMatrix(rows, cols);

    // Perform element-wise division
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            Fraction aVal = (i < r1 && j < c1) ? A[i][j] : 0;  // Treat missing elements as 0
            Fraction bVal = (i < r2 && j < c2) ? B[i][j] : 0;  // Treat missing elements as 0

            if (bVal != 0) {
                result[i][j] = aVal / bVal;
            }
            else {
                cout << "Warning: Division by zero at element [" << i + 1 << "][" << j + 1 << "]. Setting result to 0.\n";
                result[i][j] = 0;
            }
        }
    }

    return result;
}

vector<vector<Fraction>> multiplyMatrices(const vector<vector<Fraction>>& A, const vector<vector<Fraction>>& B, int r1, int c1, int r2, int c2) {
    if (c1 != r2) {
        cout << "Matrix multiplication not possible! Number of columns in A must match the number of rows in B.\n";
        return createZeroMatrix(r1, c2);
    }

    vector<vector<Fraction>> result = createZeroMatrix(r1, c2);

    for (int i = 0; i < r1; ++i) {
        for (int j = 0; j < c2; ++j) {
            for (int k = 0; k < c1; ++k) {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return result;
}

void matrixOperations() {
    int choice;
    cout << "\nSelect matrix operation:\n";
    cout << "1. Addition\n";
    cout << "2. Subtraction\n";
    cout << "3. Multiplication\n";
    cout << "4. Division (element-wise)\n";
    cout << "Enter your choice: ";
    cin >> choice;

    int r1, c1, r2, c2;
    cout << "Enter the number of rows and columns for Matrix A (e.g., 2 3): ";
    cin >> r1 >> c1;
    cout << "Enter the number of rows and columns for Matrix B (e.g., 2 3): ";
    cin >> r2 >> c2;

    vector<vector<Fraction>> A(r1, vector<Fraction>(c1)), B(r2, vector<Fraction>(c2)), result;

    inputMatrix(A, r1, c1, "A");
    inputMatrix(B, r2, c2, "B");

    switch (choice) {
    case 1:
        result = addMatrices(A, B, r1, c1, r2, c2);
        cout << "\nResult of Matrix Addition:\n";
        break;
    case 2:
        result = subtractMatrices(A, B, r1, c1, r2, c2);
        cout << "\nResult of Matrix Subtraction:\n";
        break;
    case 3:
        result = multiplyMatrices(A, B, r1, c1, r2, c2);
        cout << "\nResult of Matrix Multiplication:\n";
        break;
    case 4:
        result = divideMatrices(A, B, r1, c1, r2, c2);
        cout << "\nResult of Element-wise Matrix Division:\n";
        break;
    default:
        cout << "Invalid choice!\n";
        return;
    }

    printMatrix(result);
}

void gaussianElimination() {
    solveMatrixEquation();
}

vector<vector<Fraction>> invertMatrix(const vector<vector<Fraction>>& A, int n) {
    vector<vector<Fraction>> augmented(n, vector<Fraction>(2 * n, 0));

    // Initialize augmented matrix (A | I)
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmented[i][j] = A[i][j];
        }
        augmented[i][n + i] = 1; // Identity matrix
    }

    cout << "\nInitial Augmented Matrix (A | I):\n";
    printMatrix(augmented);

    // Perform Gaussian Elimination
    for (int i = 0; i < n; ++i) {
        // Pivoting
        Fraction maxElement = augmented[i][i].abs();
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (augmented[k][i].abs() > maxElement) {
                maxElement = augmented[k][i].abs();
                maxRow = k;
            }
        }

        if (maxRow != i) {
            swap(augmented[i], augmented[maxRow]);
            cout << "Swapped row " << i + 1 << " with row " << maxRow + 1 << " for pivoting.\n";
            printMatrix(augmented);
        }

        // Make the diagonal element 1
        Fraction diagElement = augmented[i][i];
        const Fraction epsilon = 1e-10; // Adjust as necessary
        if (diagElement.abs() < epsilon) {
            cout << "Matrix is singular and cannot be inverted.\n";
            return createZeroMatrix(n, n);
        }
        for (int j = 0; j < 2 * n; ++j) {
            augmented[i][j] /= diagElement;
        }
        cout << "Normalized row " << i + 1 << ":\n";
        printMatrix(augmented);

        // Eliminate the other elements in this column
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                Fraction factor = augmented[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
                cout << "Eliminated element in row " << k + 1 << ", column " << i + 1 << ":\n";
                printMatrix(augmented);
            }
        }
    }

    // Extract the inverse matrix from the augmented matrix
    vector<vector<Fraction>> inverse(n, vector<Fraction>(n, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverse[i][j] = augmented[i][n + j];
        }
    }

    return inverse;
}

void inverseMatrixCalculation() {
    int n;
    cout << "Enter the dimension of the square matrix (n x n): ";
    cin >> n;

    vector<vector<Fraction>> mat(n, vector<Fraction>(n));
    inputMatrix(mat, n, n, "A");

    cout << "\nOriginal Matrix A:\n";
    printMatrix(mat);

    vector<vector<Fraction>> inverse = invertMatrix(mat, n);

    // Check if inverse exists
    if (inverse.empty()) {
        cout << "Inverse matrix could not be calculated (Matrix is singular).\n";
        return;
    }

    cout << "Inverse Matrix A^(-1):\n";
    printMatrix(inverse);

    // Input B matrix
    vector<vector<Fraction>> B(n, vector<Fraction>(n));
    inputMatrix(B, n, n, "B");

    // Calculate resulting matrix X = B * A^(-1)
    vector<vector<Fraction>> X(n, vector<Fraction>(n, 0));
    X = multiplyMatrices(B, inverse, n, n, n, n);

    cout << "Resulting Matrix X = B * A^(-1):\n";
    printMatrix(X);

}

// 5. Principal Element Method (no idea)
void principalElementMethod() {
    int n;
    cout << "Enter the number of variables (n): ";
    cin >> n;

    vector<vector<Fraction>> augmented(n, vector<Fraction>(n + 1));

    cout << "\nWe are solving a system of linear equations represented as an augmented matrix (Ax = B).\n";
    cout << "Please input the coefficients of the augmented matrix (n rows, n+1 columns, where the last column is the constant vector B).\n";

    inputMatrix(augmented, n, n + 1, "Augmented Matrix");

    cout << "\nInitial Augmented Matrix:\n";
    printMatrix(augmented);

    // Principal Element Method - Forward Elimination
    for (int i = 0; i < n; ++i) {
        // Find the row with the largest element in the current column
        Fraction maxElement = augmented[i][i].abs();
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (augmented[k][i].abs() > maxElement) {
                maxElement = augmented[k][i].abs();
                maxRow = k;
            }
        }

        // Swap the current row with the row containing the principal element
        if (maxRow != i) {
            swap(augmented[i], augmented[maxRow]);
            cout << "Swapped row " << i + 1 << " with row " << maxRow + 1 << " for pivoting.\n";
            printMatrix(augmented);
        }

        // Normalize the current row by dividing by the principal element
        Fraction principalElement = augmented[i][i];
        if (principalElement == 0) {
            cout << "No unique solution exists (singular matrix).\n";
            return;
        }

        for (int j = i; j <= n; ++j) {
            augmented[i][j] /= principalElement;
        }
        cout << "Normalized row " << i + 1 << ":\n";
        printMatrix(augmented);

        // Eliminate the elements in the column below the principal element
        for (int k = i + 1; k < n; ++k) {
            Fraction factor = augmented[k][i];  // Elimination factor
            for (int j = i; j <= n; ++j) {
                augmented[k][j] -= factor * augmented[i][j];
            }
            cout << "Eliminated element in row " << k + 1 << ", column " << i + 1 << " using row " << i + 1 << ".\n";

            printMatrix(augmented);
        }
    }

    // Back Substitution to get the solution
    vector<Fraction> x(n, 0);
    for (int i = n - 1; i >= 0; --i) {
        x[i] = augmented[i][n];  // Start with the constant term
        for (int j = i + 1; j < n; ++j) {
            x[i] -= augmented[i][j] * x[j];  // Subtract the known variables
        }
    }

    cout << "Solution:\n";
    for (int i = 0; i < n; ++i) {
        cout << "x" << i + 1 << " = " << fixed << setprecision(4) << x[i] << endl;
    }
}

// Main function to display menu and call appropriate methods
int main() {
    int choice;
    do {
        cout << "\nSelect the operation (EVERYTHING DONE BY MOODLE FORMULAS EXAMPLES):\n";
        cout << "1. Matrix equation solver\n";
        cout << "2. Sum, Subtraction, Multiplication, Division\n";
        cout << "3. Gaussian Method\n";
        cout << "4. Inverse matrix calculation\n";
        cout << "5. Principal element method\n";
        cout << "0. Exit\n";
        cout << "Enter your choice: ";
        cin >> choice;

        switch (choice) {
        case 1:
            solveMatrixEquation();
            break;
        case 2:
        {
            matrixOperations();
        }
        break;
        case 3:
            gaussianElimination();
            break;
        case 4:
            inverseMatrixCalculation();
            break;
        case 5:
            principalElementMethod();
            break;
        case 0:
            cout << "Exiting...\n";
            break;
        default:
            cout << "Invalid choice! Please choose a valid option.\n";
        }
    } while (choice != 0);

    return 0;
}
