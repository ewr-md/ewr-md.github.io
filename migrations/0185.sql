ALTER TABLE chores
ADD rescheduled_next_execution_assigned_to_user_id INT;

DROP VIEW chores_current;
CREATE VIEW chores_current
AS
SELECT
	x.chore_id AS id, -- Dummy, LessQL needs an id column
	x.chore_id,
	x.chore_name,
	x.last_tracked_time,
	CASE WHEN x.rollover = 1 AND DATETIME('now', 'localtime') > x.next_estimated_execution_time THEN
		CASE WHEN IFNULL(x.track_date_only, 0) = 1 THEN
			DATETIME(STRFTIME('%Y-%m-%d', DATETIME('now', 'localtime')) || ' 23:59:59')
		ELSE
			DATETIME(STRFTIME('%Y-%m-%d', DATETIME('now', 'localtime')) || ' ' || STRFTIME('%H:%M:%S', x.next_estimated_execution_time))
		END
	ELSE
		CASE WHEN IFNULL(x.track_date_only, 0) = 1 THEN
			DATETIME(STRFTIME('%Y-%m-%d', x.next_estimated_execution_time) || ' 23:59:59')
		ELSE
			x.next_estimated_execution_time
		END
	END AS next_estimated_execution_time,
	x.track_date_only,
	x.next_execution_assigned_to_user_id,
	CASE WHEN IFNULL(x.rescheduled_date, '') != '' THEN 1 ELSE 0 END AS is_rescheduled,
	CASE WHEN IFNULL(x.rescheduled_next_execution_assigned_to_user_id, '') != '' THEN 1 ELSE 0 END AS is_reassigned
FROM (

SELECT
	h.id AS chore_id,
	h.name AS chore_name,
	MAX(l.tracked_time) AS last_tracked_time,
	CASE WHEN IFNULL(h.rescheduled_date, '') != '' THEN
		h.rescheduled_date
	ELSE
		CASE WHEN MAX(l.tracked_time) IS NULL AND h.period_type != 'manually' THEN
			h.start_date
		ELSE
			CASE h.period_type
				WHEN 'manually' THEN NULL
				WHEN 'hourly' THEN DATETIME(MAX(l.tracked_time), '+' || CAST(h.period_interval AS TEXT) || ' hour')
				WHEN 'daily' THEN DATETIME(SUBSTR(CAST(DATETIME(MAX(l.tracked_time), '+' || CAST(h.period_interval AS TEXT) || ' days') AS TEXT), 1, 11) || SUBSTR(CAST(h.start_date AS TEXT), -8))
				WHEN 'weekly' THEN (
					SELECT next
						FROM (
						SELECT 'sunday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 0') AS next
						UNION
						SELECT 'monday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 1') AS next
						UNION
						SELECT 'tuesday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 2') AS next
						UNION
						SELECT 'wednesday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 3') AS next
						UNION
						SELECT 'thursday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 4') AS next
						UNION
						SELECT 'friday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 5') AS next
						UNION
						SELECT 'saturday' AS day, DATETIME((SELECT tracked_time FROM chores_log WHERE chore_id = h.id ORDER BY tracked_time DESC LIMIT 1), '1 days', '+' || CAST((h.period_interval - 1) * 7 AS TEXT) || ' days', 'weekday 6') AS next
					)
					WHERE INSTR(period_config, day) > 0
					ORDER BY next
					LIMIT 1
				)
				WHEN 'monthly' THEN DATETIME(MAX(l.tracked_time), 'start of month', '+' || CAST(h.period_interval AS TEXT) || ' month', '+' || CAST(h.period_days - 1 AS TEXT) || ' day')
				WHEN 'yearly' THEN DATETIME(SUBSTR(CAST(DATETIME(MAX(l.tracked_time), '+' || CAST(h.period_interval AS TEXT) || ' years') AS TEXT), 1, 4) || SUBSTR(CAST(h.start_date AS TEXT), 5, 6) || SUBSTR(CAST(DATETIME(MAX(l.tracked_time), '+' || CAST(h.period_interval AS TEXT) || ' years') AS TEXT), -9))
				WHEN 'adaptive' THEN DATETIME(MAX(l.tracked_time), '+' || CAST(IFNULL((SELECT average_frequency_hours FROM chores_execution_average_frequency WHERE chore_id = h.id), 0) AS TEXT) || ' hour')
			END
		END
	END AS next_estimated_execution_time,
	h.track_date_only,
	h.rollover,
	h.next_execution_assigned_to_user_id,
	h.rescheduled_date,
	h.rescheduled_next_execution_assigned_to_user_id
FROM chores h
LEFT JOIN chores_log l
	ON h.id = l.chore_id
	AND l.undone = 0
WHERE h.active = 1
GROUP BY h.id, h.name, h.period_days
) x;
