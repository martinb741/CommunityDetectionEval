package cdTest;

import java.sql.Connection;
import java.sql.SQLException;

import org.apache.commons.dbcp2.BasicDataSource;

/**
 * This class manages database credentials and provides connection from a connection pooling system
 *
 */
public class DatabaseManager{

	private static String jdbcDriverClassName = "com.mysql.jdbc.Driver";
	private static String jdbcLogin = "recommendws";
	private static String jdbcPass = "recws!8375";
	private static String jdbcUrl = "jdbc:mysql://localhost:3306/";
	private static String jdbcSchema = "Recommender";
	private static BasicDataSource dataSource;

	public DatabaseManager() {
		dataSource = new BasicDataSource();
		dataSource.setDefaultAutoCommit(true);
		dataSource.setDriverClassName(jdbcDriverClassName);
		dataSource.setUsername(jdbcLogin);
		dataSource.setPassword(jdbcPass);
		dataSource.setUrl(jdbcUrl + jdbcSchema);
		dataSource.setValidationQuery("SELECT 1");
		dataSource.setDefaultQueryTimeout(1000);
		dataSource.setMaxConnLifetimeMillis(100000);
	}

	public Connection getConnection() throws SQLException {
		return dataSource.getConnection();
	}
	
	
	
}
