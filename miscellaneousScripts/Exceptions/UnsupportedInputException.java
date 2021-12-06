package Exceptions;

//https://rollbar.com/guides/java-throwing-exceptions/
public class UnsupportedInputException extends Exception{
    private String code;

    public UnsupportedInputException(String code, String message) {
        super(message);
        this.setCode(code);
    }

    public UnsupportedInputException(String code, String message, Throwable cause) {
        super(message, cause);
        this.setCode(code);
    }

    public String getCode() {
        return code;
    }

    public void setCode(String code) {
        this.code = code;
    }
}
